#include "global.h"
#include "build_graph.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "async_read_io.h"
#include "seq_loading_stats.h"
#include "util.h"
#include "file_util.h"

#include <pthread.h>
#include "seq_file/seq_file.h"

// Update shared_nreads in steps of 100 to reduce thread interaction
#define BUILD_GRAPH_COUNTER_STEP 100

typedef struct {
  dBGraph *db_graph;
  SeqLoadingStats *stats; // [files]
  size_t nreads;
  volatile size_t *shared_nreads;
} BuildGraphThread;

//
// Check for PCR duplicates
//

// Read start (duplicate removal during read loading)
#define db_node_has_read_start_mt(graph,node) \
        bitset_get_mt((graph)->readstrt, 2*(node).key+(node).orient)
#define db_node_set_read_start_mt(graph,node) \
        bitset_set_mt((graph)->readstrt, 2*(node).key+(node).orient)

// Returns true if start1, start2 set and reads should be added
static bool seq_reads_are_novel(read_t *r1, read_t *r2,
                                uint8_t fq_cutoff1, uint8_t fq_cutoff2,
                                uint8_t hp_cutoff, ReadMateDir matedir,
                                SeqLoadingStats *stats, dBGraph *db_graph)
{
  // Remove SAM/BAM duplicates
  if(r1->from_sam && seq_read_bam(r1)->core.flag & BAM_FDUP &&
     (r2 == NULL || (r2->from_sam && seq_read_bam(r2)->core.flag & BAM_FDUP))) {
    return false;
  }

  seq_reader_orient_mp_FF(r1, r2, matedir);

  const size_t kmer_size = db_graph->kmer_size;
  size_t start1, start2 = 0;
  bool got_kmer1 = false, got_kmer2 = false;
  BinaryKmer bkmer1, bkmer2;
  dBNode node1 = DB_NODE_INIT, node2 = DB_NODE_INIT;

  start1 = seq_contig_start(r1, 0, kmer_size, fq_cutoff1, hp_cutoff);
  got_kmer1 = (start1 < r1->seq.end);

  if(r2) {
    start2 = seq_contig_start(r2, 0, kmer_size, fq_cutoff2, hp_cutoff);
    got_kmer2 = (start2 < r2->seq.end);
  }

  bool found1 = false, found2 = false;

  // Look up first kmer
  if(got_kmer1) {
    bkmer1 = binary_kmer_from_str(r1->seq.b + start1, kmer_size);
    node1 = db_graph_find_or_add_node_mt(db_graph, bkmer1, &found1);
  }

  // Look up second kmer
  if(got_kmer2) {
    bkmer2 = binary_kmer_from_str(r2->seq.b + start2, kmer_size);
    node2 = db_graph_find_or_add_node_mt(db_graph, bkmer2, &found2);
  }

  size_t num_kmers_novel = !found1 + !found2;
  __sync_fetch_and_add((volatile size_t*)&stats->num_kmers_novel, num_kmers_novel);

  // Each read gives no kmer or a duplicate kmer
  // used find_or_insert so if we have a kmer we have a graph node
  if((!got_kmer1 || db_node_has_read_start_mt(db_graph, node1)) &&
     (!got_kmer2 || db_node_has_read_start_mt(db_graph, node2)))
  {
    return false;
  }

  // Read is novel
  if(got_kmer1) (void)db_node_set_read_start_mt(db_graph, node1);
  if(got_kmer2) (void)db_node_set_read_start_mt(db_graph, node2);

  return true;
}


//
// Add to the de bruijn graph
//

static inline dBNode _find_or_insert(dBGraph *db_graph,
                                     BinaryKmer bkmer, size_t colour,
                                     bool must_exist_in_graph, bool *found)
{
  dBNode node;
  if(must_exist_in_graph)
  {
    // Doesn't have to be threadsafe find_mt, since we are not adding
    node = db_graph_find_node(db_graph, bkmer);
    *found = (node.key != HASH_NOT_FOUND);
    if(*found) db_graph_update_node_mt(db_graph, node, colour);
  }
  else
  {
    node = db_graph_find_or_add_node_mt(db_graph, bkmer, found);
    db_graph_update_node_mt(db_graph, node, colour);
  }
  return node;
}

// Threadsafe
// Sequence must be entirely ACGT and len >= kmer_size
// Returns number of non-novel kmers seen
size_t build_graph_from_str_mt(dBGraph *db_graph, size_t colour,
                               const char *seq, size_t len,
                               bool must_exist_in_graph)
{
  ctx_assert(len >= db_graph->kmer_size);
  const size_t kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer;
  Nucleotide nuc;
  dBNode prev, curr;
  size_t i, num_nonnovel_kmers = 0;
  size_t edge_col = db_graph->num_edge_cols == 1 ? 0 : colour;
  bool found;

  bkmer = binary_kmer_from_str(seq, kmer_size);
  prev = _find_or_insert(db_graph, bkmer, colour, must_exist_in_graph, &found);
  num_nonnovel_kmers += found;

  for(i = kmer_size; i < len; i++, prev = curr)
  {
    nuc = dna_char_to_nuc(seq[i]);
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    curr = _find_or_insert(db_graph, bkmer, colour, must_exist_in_graph, &found);
    if(prev.key != HASH_NOT_FOUND && curr.key != HASH_NOT_FOUND)
      db_graph_add_edge_mt(db_graph, edge_col, prev, curr);
    num_nonnovel_kmers += found;
  }

  return num_nonnovel_kmers;
}

// Already found a start position
// Stats must be private to this thread
static void load_read(const read_t *r, uint8_t qual_cutoff, uint8_t hp_cutoff,
                      bool must_exist_in_graph, Colour colour,
                      SeqLoadingStats *stats, dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  size_t contig_start, contig_end, contig_len;
  size_t num_contigs = 0, search_start = 0, num_nonnovel_kmers;

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qual_cutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qual_cutoff, hp_cutoff, &search_start);

    contig_len = contig_end - contig_start;
    num_nonnovel_kmers = build_graph_from_str_mt(db_graph, colour,
                                                 r->seq.b+contig_start, contig_len,
                                                 must_exist_in_graph);

    size_t contig_kmers = contig_len + 1 - kmer_size;
    size_t num_novel_kmers = contig_kmers - num_nonnovel_kmers;

    stats->total_bases_loaded += contig_len;
    if(must_exist_in_graph) {
      stats->num_kmers_loaded += num_nonnovel_kmers;
    } else {
      stats->num_kmers_loaded += contig_kmers;
      stats->num_kmers_novel += num_novel_kmers;
    }
    num_contigs++;
  }

  stats->contigs_parsed += num_contigs;
  stats->num_good_reads += (num_contigs > 0);
  stats->num_bad_reads += (num_contigs == 0);
}

// Stats must be private to this thread
void build_graph_from_reads_mt(read_t *r1, read_t *r2,
                               uint8_t fq_offset1, uint8_t fq_offset2,
                               const SeqLoadingPrefs *prefs,
                               SeqLoadingStats *stats,
                               dBGraph *db_graph)
{
  ctx_assert(!prefs->must_exist_in_graph || !prefs->remove_pcr_dups);
  // status("r1: '%s' '%s'", r1->name.b, r1->seq.b);
  // if(r2) status("r2: '%s' '%s'", r2->name.b, r2->seq.b);

  uint8_t fq_cutoff1 = prefs->fq_cutoff, fq_cutoff2 = prefs->fq_cutoff;

  if(prefs->fq_cutoff) {
    fq_cutoff1 += fq_offset1;
    fq_cutoff2 += fq_offset2;
  }

  size_t total_bases = r1->seq.end + (r2 ? r2->seq.end : 0);
  stats->total_bases_read += total_bases;

  if(r2) stats->num_pe_reads += 2;
  else   stats->num_se_reads += 1;

  // printf(">%s %zu\n", r1->name.b, colour);

  if(prefs->remove_pcr_dups && !seq_reads_are_novel(r1, r2,
                                                   fq_cutoff1, fq_cutoff2,
                                                   prefs->hp_cutoff, prefs->matedir,
                                                   stats, db_graph))
  {
    if(r2) stats->num_dup_pe_pairs++;
    else   stats->num_dup_se_reads++;
  }
  else {
    load_read(r1, fq_cutoff1, prefs->hp_cutoff, prefs->must_exist_in_graph,
              prefs->colour, stats, db_graph);
    if(r2) load_read(r2, fq_cutoff2, prefs->hp_cutoff, prefs->must_exist_in_graph,
                     prefs->colour, stats, db_graph);
  }
}

static void add_reads_to_graph(AsyncIOData *data, size_t threadid, void *ptr)
{
  (void)threadid;
  BuildGraphThread *wrkr = (BuildGraphThread*)ptr;
  const BuildGraphTask *task = (BuildGraphTask*)data->ptr;
  read_t *r2 = data->r2.name.end == 0 && data->r2.seq.end == 0 ? NULL : &data->r2;

  build_graph_from_reads_mt(&data->r1, r2,
                            data->fq_offset1, data->fq_offset2,
                            &task->prefs, wrkr->stats,
                            wrkr->db_graph);

  // Print progress
  wrkr->nreads++;
  if(wrkr->nreads >= BUILD_GRAPH_COUNTER_STEP) {
    // Update shared counter
    size_t n = __sync_fetch_and_add(wrkr->shared_nreads, wrkr->nreads);
    // if n .. n+wrkr->nreads
    ctx_update2("BuildGraph", n, n+wrkr->nreads, CTX_UPDATE_REPORT_RATE);
    wrkr->nreads = 0;
  }
}

// One thread used per input file, nthreads used to add reads to graph
void build_graph(dBGraph *db_graph, BuildGraphTask *files,
                 size_t nfiles, size_t nthreads)
{
  ctx_assert(db_graph->bktlocks != NULL);

  // Start async io reading
  AsyncIOInput *async_tasks = ctx_malloc(nfiles * sizeof(AsyncIOInput));
  size_t i, f;

  for(f = 0; f < nfiles; f++) {
    files[f].idx = f;
    files[f].files.ptr = &files[f];
    memcpy(&async_tasks[f], &files[f].files, sizeof(AsyncIOInput));
  }

  BuildGraphThread *threads = ctx_calloc(nthreads, sizeof(BuildGraphThread));
  size_t total_nreads = 0;

  for(i = 0; i < nthreads; i++) {
    threads[i].stats = ctx_calloc(nfiles, sizeof(SeqLoadingStats));
    threads[i].db_graph = db_graph;
    threads[i].shared_nreads = &total_nreads;
  }

  asyncio_run_pool(async_tasks, nfiles, add_reads_to_graph,
                   threads, nthreads, sizeof(BuildGraphThread));

  // Merge stats
  for(i = 0; i < nthreads; i++) {
    for(f = 0; f < nfiles; f++)
      seq_loading_stats_merge(&files[f].stats, &threads[i].stats[f]);
    ctx_free(threads[i].stats);
  }
  ctx_free(threads);
  ctx_free(async_tasks);

  // Copy stats into ginfo
  size_t max_col = 0;
  for(f = 0; f < nfiles; f++) {
    max_col = MAX2(max_col, files[f].prefs.colour);
    graph_info_update_stats(&db_graph->ginfo[files[f].prefs.colour], &files[f].stats);
  }

  db_graph->num_of_cols_used = MAX2(db_graph->num_of_cols_used, max_col+1);
}

// One thread used per input file, nthreads used to add reads to graph
// Updates ginfo
void build_graph_from_seq(dBGraph *db_graph,
                          seq_file_t **files, size_t nfiles,
                          size_t nthreads,
                          size_t colour)
{
  size_t i;
  BuildGraphTask *tasks = ctx_calloc(nfiles, sizeof(BuildGraphTask));

  for(i = 0; i < nfiles; i++) {
    AsyncIOInput input = {.file1 = files[i], .file2 = NULL,
                          .fq_offset = 0, .interleaved = false};

    BuildGraphTask tmp = {.files = input,
                          .prefs = SEQ_LOADING_PREFS_INIT,
                          .stats = SEQ_LOADING_STATS_INIT,
                          .idx = 0};

    tmp.prefs.colour = colour;
    seq_loading_stats_init(&tmp.stats);
    memcpy(&tasks[i], &tmp, sizeof(tmp));
  }

  build_graph(db_graph, tasks, nfiles, nthreads);

  ctx_free(tasks);
}

void build_graph_task_print(const BuildGraphTask *task)
{
  const AsyncIOInput *io = &task->files;
  const SeqLoadingPrefs *prefs = &task->prefs;

  char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

  if(io->fq_offset > 0) sprintf(fqOffset, "%u", io->fq_offset);
  if(prefs->fq_cutoff > 0) sprintf(fqCutoff, "%u", prefs->fq_cutoff);
  if(prefs->hp_cutoff > 0) sprintf(hpCutoff, "%u", prefs->hp_cutoff);

  status("[task] %s%s%s; FASTQ offset: %s, threshold: %s; "
         "cut homopolymers: %s; remove PCR duplicates: %s; colour: %zu\n",
         futil_inpath_str(io->file1->path),
         io->file2 ? ", " : "",
         io->file2 ? futil_inpath_str(io->file2->path) : "",
         fqOffset, fqCutoff, hpCutoff, prefs->remove_pcr_dups ? "yes" : "no",
         prefs->colour);
}

void build_graph_task_print_stats(const BuildGraphTask *task)
{
  const SeqLoadingStats *stats = &task->stats;
  const SeqLoadingPrefs *prefs = &task->prefs;
  const AsyncIOInput *io = &task->files;

  status("[task] input: %s%s%s colour: %zu",
         futil_inpath_str(io->file1->path), io->file2 ? ", " : "",
         io->file2 ? futil_inpath_str(io->file2->path) : "", prefs->colour);

  char se_reads_str[50], pe_reads_str[50];
  char good_reads_str[50], bad_reads_str[50];
  char dup_se_reads_str[50], dup_pe_pairs_str[50];
  char bases_read_str[50], bases_loaded_str[50];
  char num_contigs_str[50], num_kmers_loaded_str[50], num_kmers_novel_str[50];

  ulong_to_str(stats->num_se_reads, se_reads_str);
  ulong_to_str(stats->num_pe_reads, pe_reads_str);
  ulong_to_str(stats->num_good_reads, good_reads_str);
  ulong_to_str(stats->num_bad_reads, bad_reads_str);
  ulong_to_str(stats->num_dup_se_reads, dup_se_reads_str);
  ulong_to_str(stats->num_dup_pe_pairs, dup_pe_pairs_str);
  ulong_to_str(stats->total_bases_read, bases_read_str);
  ulong_to_str(stats->total_bases_loaded, bases_loaded_str);
  ulong_to_str(stats->contigs_parsed, num_contigs_str);
  ulong_to_str(stats->num_kmers_loaded, num_kmers_loaded_str);
  ulong_to_str(stats->num_kmers_novel, num_kmers_novel_str);

  status("  SE reads: %s  PE reads: %s", se_reads_str, pe_reads_str);
  status("  good reads: %s  bad reads: %s", good_reads_str, bad_reads_str);
  status("  dup SE reads: %s  dup PE pairs: %s", dup_se_reads_str, dup_pe_pairs_str);
  status("  bases read: %s  bases loaded: %s", bases_read_str, bases_loaded_str);
  status("  num contigs: %s  num kmers: %s novel kmers: %s",
         num_contigs_str, num_kmers_loaded_str, num_kmers_novel_str);
}
