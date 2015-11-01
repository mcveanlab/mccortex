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

typedef struct
{
  dBGraph *const db_graph;
  size_t rcounter; // counter of entries taken from the pool
} BuildGraphData;

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

// Threadsafe
// Sequence must be entirely ACGT and len >= kmer_size
// Returns number of novel kmers loaded
size_t build_graph_from_str_mt(dBGraph *db_graph, size_t colour,
                               const char *seq, size_t len)
{
  ctx_assert(len >= db_graph->kmer_size);
  const size_t kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer;
  Nucleotide nuc;
  dBNode prev, curr;
  size_t i, num_novel_kmers = 0;
  size_t edge_col = db_graph->num_edge_cols == 1 ? 0 : colour;
  bool found;

  bkmer = binary_kmer_from_str(seq, kmer_size);
  prev = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
  db_graph_update_node_mt(db_graph, prev, colour);
  num_novel_kmers += !found;

  for(i = kmer_size; i < len; i++)
  {
    nuc = dna_char_to_nuc(seq[i]);
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    curr = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
    db_graph_update_node_mt(db_graph, curr, colour);
    db_graph_add_edge_mt(db_graph, edge_col, prev, curr);
    num_novel_kmers += !found;
    prev = curr;
  }

  return num_novel_kmers;
}

// Already found a start position
static void load_read(const read_t *r, uint8_t qual_cutoff, uint8_t hp_cutoff,
                      SeqLoadingStats *stats, Colour colour, dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  size_t contig_start, contig_end, contig_len;
  size_t num_contigs = 0, search_start = 0, num_novel_kmers;

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qual_cutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qual_cutoff, hp_cutoff, &search_start);

    contig_len = contig_end - contig_start;
    num_novel_kmers = build_graph_from_str_mt(db_graph, colour,
                                              r->seq.b+contig_start, contig_len);

    size_t contig_kmers = contig_len + 1 - kmer_size;
    __sync_fetch_and_add((volatile size_t*)&stats->total_bases_loaded, contig_len);
    __sync_fetch_and_add((volatile size_t*)&stats->num_kmers_loaded, contig_kmers);
    __sync_fetch_and_add((volatile size_t*)&stats->num_kmers_novel, num_novel_kmers);
    num_contigs++;
  }

  __sync_fetch_and_add((volatile size_t*)&stats->contigs_parsed, num_contigs);
  __sync_fetch_and_add((volatile size_t*)&stats->num_good_reads, num_contigs > 0);
  __sync_fetch_and_add((volatile size_t*)&stats->num_bad_reads, num_contigs == 0);
}

void build_graph_from_reads_mt(read_t *r1, read_t *r2,
                               uint8_t fq_offset1, uint8_t fq_offset2,
                               uint8_t fq_cutoff, uint8_t hp_cutoff,
                               bool remove_pcr_dups, ReadMateDir matedir,
                               SeqLoadingStats *stats, size_t colour,
                               dBGraph *db_graph)
{
  // status("r1: '%s' '%s'", r1->name.b, r1->seq.b);
  // if(r2) status("r2: '%s' '%s'", r2->name.b, r2->seq.b);

  uint8_t fq_cutoff1 = fq_cutoff, fq_cutoff2 = fq_cutoff;

  if(fq_cutoff) {
    fq_cutoff1 += fq_offset1;
    fq_cutoff2 += fq_offset2;
  }

  size_t total_bases = r1->seq.end + (r2 ? r2->seq.end : 0);
  __sync_fetch_and_add((volatile size_t*)&stats->total_bases_read, total_bases);

  if(r2) __sync_fetch_and_add((volatile size_t*)&stats->num_pe_reads, 2);
  else   __sync_fetch_and_add((volatile size_t*)&stats->num_se_reads, 1);

  // printf(">%s %zu\n", r1->name.b, colour);

  if(remove_pcr_dups && !seq_reads_are_novel(r1, r2,
                                             fq_cutoff1, fq_cutoff2, hp_cutoff,
                                             matedir, stats, db_graph))
  {
    if(r2) __sync_fetch_and_add((volatile size_t*)&stats->num_dup_pe_pairs, 1);
    else   __sync_fetch_and_add((volatile size_t*)&stats->num_dup_se_reads, 1);
  }
  else {
    load_read(r1, fq_cutoff1, hp_cutoff, stats, colour, db_graph);
    if(r2) load_read(r2, fq_cutoff2, hp_cutoff, stats, colour, db_graph);
  }
}

static void add_reads_to_graph(AsyncIOData *data, void *ptr)
{
  BuildGraphData *wrkr = (BuildGraphData*)ptr;
  BuildGraphTask *task = (BuildGraphTask*)data->ptr;
  read_t *r2 = data->r2.name.end == 0 && data->r2.seq.end == 0 ? NULL : &data->r2;

  build_graph_from_reads_mt(&data->r1, r2,
                            data->fq_offset1, data->fq_offset2,
                            task->fq_cutoff, task->hp_cutoff,
                            task->remove_pcr_dups, task->matedir,
                            &task->stats,
                            task->colour, wrkr->db_graph);

  // Print progress
  size_t n = __sync_add_and_fetch((volatile size_t*)&wrkr->rcounter, 1);
  ctx_update("BuildGraph", n);
}

// One thread used per input file, num_build_threads used to add reads to graph
void build_graph(dBGraph *db_graph, BuildGraphTask *files,
                 size_t num_files, size_t num_build_threads)
{
  ctx_assert(db_graph->bktlocks != NULL);

  // Start async io reading
  AsyncIOInput *async_tasks = ctx_malloc(num_files * sizeof(AsyncIOInput));
  size_t f;

  for(f = 0; f < num_files; f++) {
    files[f].idx = f;
    files[f].files.ptr = &files[f];
    memcpy(&async_tasks[f], &files[f].files, sizeof(AsyncIOInput));
  }

  BuildGraphData global_mem = {.db_graph = db_graph, .rcounter = 0};

  asyncio_run_pool(async_tasks, num_files, add_reads_to_graph,
                   &global_mem, num_build_threads, 0);

  ctx_free(async_tasks);

  // Copy stats into ginfo
  size_t max_col = 0;
  for(f = 0; f < num_files; f++) {
    max_col = MAX2(max_col, files[f].colour);
    graph_info_update_stats(&db_graph->ginfo[files[f].colour], &files[f].stats);
  }

  db_graph->num_of_cols_used = MAX2(db_graph->num_of_cols_used, max_col+1);
}

// One thread used per input file, num_build_threads used to add reads to graph
// Updates ginfo
void build_graph_from_seq(dBGraph *db_graph,
                          seq_file_t **files, size_t num_files,
                          size_t num_build_threads,
                          size_t colour)
{
  size_t i;
  BuildGraphTask *tasks = ctx_calloc(num_files, sizeof(BuildGraphTask));

  for(i = 0; i < num_files; i++) {
    AsyncIOInput input = {.file1 = files[i], .file2 = NULL,
                          .fq_offset = 0, .interleaved = false};

    BuildGraphTask tmp = {.files = input, .fq_cutoff = 0, .hp_cutoff = 0,
                          .matedir = READPAIR_FR, .colour = colour,
                          .remove_pcr_dups = false};

    seq_loading_stats_init(&tmp.stats);
    memcpy(&tasks[i], &tmp, sizeof(tmp));
  }

  build_graph(db_graph, tasks, num_files, num_build_threads);

  ctx_free(tasks);
}

void build_graph_task_print(const BuildGraphTask *task)
{
  const AsyncIOInput *io = &task->files;
  char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

  if(io->fq_offset > 0) sprintf(fqOffset, "%u", io->fq_offset);
  if(task->fq_cutoff > 0) sprintf(fqCutoff, "%u", task->fq_cutoff);
  if(task->hp_cutoff > 0) sprintf(hpCutoff, "%u", task->hp_cutoff);

  status("[task] %s%s%s; FASTQ offset: %s, threshold: %s; "
         "cut homopolymers: %s; remove PCR duplicates: %s; colour: %zu\n",
         futil_inpath_str(io->file1->path),
         io->file2 ? ", " : "",
         io->file2 ? futil_inpath_str(io->file2->path) : "",
         fqOffset, fqCutoff, hpCutoff, task->remove_pcr_dups ? "yes" : "no",
         task->colour);
}

void build_graph_task_print_stats(const BuildGraphTask *task)
{
  const SeqLoadingStats stats = task->stats;
  const AsyncIOInput *io = &task->files;

  status("[task] input: %s%s%s colour: %zu",
         futil_inpath_str(io->file1->path), io->file2 ? ", " : "",
         io->file2 ? futil_inpath_str(io->file2->path) : "", task->colour);

  char se_reads_str[50], pe_reads_str[50];
  char good_reads_str[50], bad_reads_str[50];
  char dup_se_reads_str[50], dup_pe_pairs_str[50];
  char bases_read_str[50], bases_loaded_str[50];
  char num_contigs_str[50], num_kmers_loaded_str[50], num_kmers_novel_str[50];

  ulong_to_str(stats.num_se_reads, se_reads_str);
  ulong_to_str(stats.num_pe_reads, pe_reads_str);
  ulong_to_str(stats.num_good_reads, good_reads_str);
  ulong_to_str(stats.num_bad_reads, bad_reads_str);
  ulong_to_str(stats.num_dup_se_reads, dup_se_reads_str);
  ulong_to_str(stats.num_dup_pe_pairs, dup_pe_pairs_str);
  ulong_to_str(stats.total_bases_read, bases_read_str);
  ulong_to_str(stats.total_bases_loaded, bases_loaded_str);
  ulong_to_str(stats.contigs_parsed, num_contigs_str);
  ulong_to_str(stats.num_kmers_loaded, num_kmers_loaded_str);
  ulong_to_str(stats.num_kmers_novel, num_kmers_novel_str);

  status("  SE reads: %s  PE reads: %s", se_reads_str, pe_reads_str);
  status("  good reads: %s  bad reads: %s", good_reads_str, bad_reads_str);
  status("  dup SE reads: %s  dup PE pairs: %s", dup_se_reads_str, dup_pe_pairs_str);
  status("  bases read: %s  bases loaded: %s", bases_read_str, bases_loaded_str);
  status("  num contigs: %s  num kmers: %s novel kmers: %s",
         num_contigs_str, num_kmers_loaded_str, num_kmers_novel_str);
}
