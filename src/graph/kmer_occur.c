#include "global.h"
#include "kmer_occur.h"
#include "seq_reader.h"
#include "util.h"
#include "db_node.h"

//
// This file provides a datastore for loading sequences and recording where
// each kmer occurs in the sequences. Used in breakpoint_caller.c.
//

// Mostly used for debugging
void korun_print(KOccurRun run, size_t kmer_size, FILE *fout)
{
  const char strand[2] = "+-";
  size_t start, end;
  if(run.strand == STRAND_PLUS) {
    start = run.first;
    end = run.last+kmer_size-1;
  } else {
    start = run.first+kmer_size-1;
    end = run.last;
  }
  fprintf(fout, "offset:%zu:chromid:%zu:%zu-%zu:%c",
          (size_t)run.qoffset, (size_t)run.chrom,
          start, end, strand[run.strand]);
}

// Get string representation of multiple runs, comma separated
//   e.g. "chromid:1:17-5:-, chromid:1:37-47:+"
// Does not print new line
// Mostly used for debugging
void koruns_print(KOccurRun *runs, size_t n, size_t kmer_size, FILE *fout)
{
  size_t i;
  if(n == 0) return;
  korun_print(runs[0], kmer_size, fout);
  for(i = 1; i < n; i++) {
    printf(", ");
    korun_print(runs[i], kmer_size, fout);
  }
}

// Sort KOccurRun by qoffset, chrom, first, last then strand
static int korun_cmp_qoffset(const void *aa, const void *bb)
{
  const KOccurRun *a = (const KOccurRun*)aa, *b = (const KOccurRun*)bb;
  if(a->qoffset != b->qoffset) return a->qoffset < b->qoffset ? -1 : 1;
  if(a->chrom != b->chrom) return a->chrom < b->chrom ? -1 : 1;
  if(a->first != b->first) return a->first < b->first ? -1 : 1;
  if(a->last != b->last) return a->last < b->last ? -1 : 1;
  return 0;
}

void koruns_sort_by_qoffset(KOccurRun *runs, size_t n)
{
  qsort(runs, n, sizeof(KOccurRun), korun_cmp_qoffset);
}

// Malloc's and returns chrom list - remember to free
// Returns NULL if num_reads == 0
static void generate_chrom_list(KOGraph *kograph,
                                const read_t *reads, size_t num_reads)
{
  kograph->chrom_name_buf = NULL;
  kograph->nchroms = num_reads;
  kograph->chroms = NULL;

  if(num_reads == 0) return;

  kograph->chroms = ctx_malloc(num_reads*sizeof(KOChrom));
  size_t i, name_len_sum = 0;

  // Concatenate names into one string
  for(i = 0; i < num_reads; i++) name_len_sum += reads[i].name.end + 1;
  kograph->chrom_name_buf = ctx_malloc(name_len_sum);
  char *names_concat = kograph->chrom_name_buf;

  for(i = 0; i < num_reads; i++)
  {
    kograph->chroms[i] = (KOChrom){.id = i,
                                   .length = reads[i].seq.end,
                                   .name = names_concat};

    memcpy(names_concat, reads[i].name.b, reads[i].name.end);
    names_concat[reads[i].name.end] = '\0';
    names_concat += reads[i].name.end+1;
  }
}

// Add missing kmers and edges to the graph whilst keeping track of the count
// of how many times each kmer in the graph is seen in sequence (with klists)
// Returns number of kmers added to the graph
// Threadsafe
static inline size_t add_ref_seq_to_graph_mt(const char *seq, size_t len,
                                             KONodeList *klists,
                                             dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  size_t i, num_novel_kmers = 0;
  BinaryKmer bkmer;
  Nucleotide nuc;
  dBNode prev, curr;
  bool found;

  ctx_assert(len >= kmer_size);

  bkmer = binary_kmer_from_str(seq, kmer_size);
  prev = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
  __sync_fetch_and_add((volatile uint32_t*)&klists[prev.key].count, 1); // count++
  num_novel_kmers += !found;

  for(i = kmer_size; i < len; i++, prev = curr)
  {
    nuc = dna_char_to_nuc(seq[i]);
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    curr = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
    __sync_fetch_and_add((volatile uint32_t*)&klists[curr.key].count, 1); // count++
    db_graph_add_edge_mt(db_graph, 0, prev, curr);
    num_novel_kmers += !found;
  }

  return num_novel_kmers;
}

// Add missing kmers and edges to the graph whilst keeping track of the count
// of how many times each kmer in the graph is seen in sequence (with klists)
// Returns number of kmers added to the graph
// Threadsafe
static inline size_t add_ref_read_to_graph_mt(const read_t *r,
                                              KONodeList *klists,
                                              dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  size_t contig_start, contig_end, contig_len;
  size_t search_start = 0, num_novel_kmers = 0;

  if(r->seq.end < kmer_size) return 0;

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         0, 0)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size, 0, 0, &search_start);

    contig_len = contig_end - contig_start;
    num_novel_kmers += add_ref_seq_to_graph_mt(r->seq.b+contig_start, contig_len,
                                               klists, db_graph);
  }

  return num_novel_kmers;
}

// Same as above but don't add missing kmers
// Threadsafe
static inline void bkmer_update_counts_find_mt(BinaryKmer bkmer,
                                               KONodeList *klists,
                                               const dBGraph *db_graph)
{
  hkey_t hkey = hash_table_find(&db_graph->ht, bkmer);
  if(hkey != HASH_NOT_FOUND)
    __sync_fetch_and_add((volatile uint32_t*)&klists[hkey].count, 1); // count++
}

struct ReadUpdateCounts {
  const read_t *r;
  KONodeList *klists;
  bool add_missing_kmers;
  dBGraph *db_graph;
};

// Multithreaded core function to store kmer occurances
// Also adds read to the graph if `add_missing_kmers` is true
static void read_update_counts(void *arg)
{
  struct ReadUpdateCounts data = *(struct ReadUpdateCounts*)arg;
  const read_t *r = data.r;

  if(data.add_missing_kmers) {
    add_ref_read_to_graph_mt(r, data.klists, data.db_graph);
  }
  else {
    LoadingStats stats = LOAD_STATS_INIT_MACRO;
    READ_TO_BKMERS(r, data.db_graph->kmer_size, 0, 0, &stats,
                   bkmer_update_counts_find_mt, data.klists, data.db_graph);
  }
}

static void bkmer_store_kmer_pos(BinaryKmer bkmer, const dBGraph *db_graph,
                                 KONodeList *klists, KOccur *koccurs,
                                 size_t chrom_id, uint64_t offset)
{
  // bkmers were already added to graph -> don't need to find_or_insert
  // if missing kmers weren't added then kmer might be missing -> skip
  dBNode node = db_graph_find(db_graph, bkmer);

  if(node.key != HASH_NOT_FOUND)
  {
    KONodeList *kl = &klists[node.key];
    koccurs[kl->start+kl->count] = (KOccur){.chrom = chrom_id,
                                            .offset = offset,
                                            .orient = node.orient};
    kl->count++;
  }
}

static void read_store_kmer_pos(const read_t *r, size_t chrom_id,
                                KONodeList *klists, KOccur *koccurs,
                                const dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  LoadingStats stats = LOAD_STATS_INIT_MACRO;
  READ_TO_BKMERS(r, kmer_size, 0, 0, &stats, bkmer_store_kmer_pos,
                 db_graph, klists, koccurs, chrom_id, _offset);
}

static void load_reads_count_kmers(const read_t *reads, size_t num_reads,
                                   bool add_missing_kmers, size_t num_threads,
                                   KONodeList *klists,
                                   dBGraph *db_graph)
{
  if(!num_reads) return;

  // 1. Loop through reads, add to graph and record kmer counts
  struct ReadUpdateCounts *updates;
  size_t i;

  updates = ctx_malloc(num_reads * sizeof(struct ReadUpdateCounts));

  for(i = 0; i < num_reads; i++) {
    updates[i] = (struct ReadUpdateCounts){.r = &reads[i],
                                           .klists = klists,
                                           .add_missing_kmers = add_missing_kmers,
                                           .db_graph = db_graph};
  }

  util_run_threads(updates, num_reads, sizeof(struct ReadUpdateCounts),
                   num_threads, read_update_counts);

  ctx_free(updates);
}

// BEWARE: We add the reads to the graph if add_missing_kmers is true
// db_graph->col_edges can be NULL even if we are adding kmers
KOGraph kograph_create(const read_t *reads, size_t num_reads,
                       bool add_missing_kmers, size_t num_threads,
                       dBGraph *db_graph)
{
  size_t i;

  status("Adding reference annotations to the graph using %zu thread%s",
         num_threads, util_plural_str(num_threads));

  // If we are adding nodes, only have edges in one colour
  //  - it gets confusing otherwise (which colour would we add edges to?)
  ctx_assert(!add_missing_kmers || db_graph->num_edge_cols <= 1);
  ctx_assert(!add_missing_kmers || db_graph->bktlocks != NULL);
  ctx_assert(sizeof(KONodeList) == 12);

  // Check number of reads doesn't exceed max limit
  if(num_reads > KMER_OCCUR_MAX_CHROMS)
    die("More chromosomes than permitted (%zu > %i)",
        num_reads, KMER_OCCUR_MAX_CHROMS);

  // Check no read is too long
  for(i = 0; i < num_reads; i++)
    if(reads[i].seq.end > KMER_OCCUR_MAX_LEN)
      die("Read longer than limit (%zu > %zu; %zu: '%s')",
          reads[i].seq.end, KMER_OCCUR_MAX_LEN, i, reads[i].name.b);

  KOGraph kograph;
  memset(&kograph, 0, sizeof(KOGraph)); // initialise

  generate_chrom_list(&kograph, reads, num_reads);

  kograph.klists = ctx_calloc(db_graph->ht.capacity, sizeof(KONodeList));

  // 1. Loop through reads, add to graph and record kmer counts
  load_reads_count_kmers(reads, num_reads, add_missing_kmers, num_threads,
                         kograph.klists, db_graph);

  // 2. allocate a list for each kmer (some of length zero)
  uint64_t offset = 0, total_read_length = 0;

  for(i = 0; i < db_graph->ht.capacity; i++)
  {
    kograph.klists[i].start = offset;
    offset += kograph.klists[i].count;
    kograph.klists[i].count = 0;
  }

  // Sum lengths of reads
  for(i = 0; i < num_reads; i++)
    total_read_length += reads[i].seq.end;

  // offset should be less than sum of read lengths, which is why we wait to
  // parse all reads before mallocing to reduce memory footprint
  ctx_assert(total_read_length == 0 || offset < total_read_length);

  if(offset == 0) kograph.koccurs = NULL;
  else {
    kograph.koccurs = ctx_malloc(offset * sizeof(KOccur));

    // 3. Loop through reads, record kmer pos
    for(i = 0; i < num_reads; i++) {
      read_store_kmer_pos(&reads[i], i, kograph.klists, kograph.koccurs, db_graph);
    }
  }

  return kograph;
}

void kograph_free(KOGraph kograph)
{
  ctx_free(kograph.chrom_name_buf);
  ctx_free(kograph.chroms);
  ctx_free(kograph.klists);
  ctx_free(kograph.koccurs);
}


//
// Check for a run of kmers in the reference genome
//

// `run` should be of length num0+num1
// `pickup` if true, create new paths for kmers not used in a path
// `qoffset` is only used is pickup is true, and is the offset in the query
// Returns length of run
static size_t korun_extend(KOccurRun *kolist0, size_t num0,
                           dBNode node,
                           const KOccur *kolist1, size_t num1,
                           KOccurRun *korun, bool pickup, size_t qoffset)
{
  size_t i, j, k = 0, starti = 0, next;
  bool used, strand;

  // By looping over kolist1 we return a sorted list (by chrom and offset)
  // takes ~ num0 + num1*3 time
  for(j = 0; j < num1; j++)
  {
    while(starti < num0 &&
          (kolist0[starti].chrom < kolist1[j].chrom ||
           kolist0[starti].last+1 < kolist1[j].offset))
    {
      starti++;
    }

    used = false;
    strand = kolist1[j].orient == node.orient ? STRAND_PLUS : STRAND_MINUS;

    for(i = starti; i < num0; i++)
    {
      if(kolist0[i].chrom > kolist1[j].chrom ||
         kolist0[i].last > kolist1[j].offset+1)
      {
        break;
      }
      else
      {
        // Sanity checks
        ctx_assert(kolist0[i].chrom == kolist1[j].chrom);
        ctx_assert(kolist0[i].last+1 == kolist1[j].offset ||
                   kolist0[i].last == kolist1[j].offset+1);

        next = strand == STRAND_PLUS ? kolist0[i].last+1 : kolist0[i].last-1;

        if(next == kolist1[j].offset) {
          korun[k++] = (KOccurRun){.chrom = kolist0[i].chrom,
                                   .first = kolist0[i].first,
                                   .last = kolist1[j].offset,
                                   .qoffset = kolist0[i].qoffset,
                                   .strand = strand,
                                   .used = false};

          kolist0[i].used = true;
          used = true;
        }
      }
    }

    if(pickup && !used)
    {
      // Start new kmer run
      korun[k++] = (KOccurRun){.chrom = kolist1[j].chrom,
                               .first = kolist1[j].offset,
                               .last = kolist1[j].offset,
                               .qoffset = qoffset,
                               .strand = strand,
                               .used = false};
    }
  }

  return k;
}

// Filter regions down to only those that stretch the whole distance
// Does not reset either korun or runs_ended - only adds to runs_ended
// korun can only get shorter as KOccurRuns finish
// `qoffset` is used for offset of new runs starting at nodes[0]
void kograph_filter_extend(KOGraph kograph,
                           const dBNode *nodes, size_t num_nodes, bool forward,
                           size_t min_len, size_t qoffset,
                           KOccurRunBuffer *korun,
                           KOccurRunBuffer *runs_ended,
                           bool pickup_at_first_node)
{
  const KOccur *kolist;
  KOccurRun *runs0, *runs1;
  size_t i, j, num_occur, nruns0, nruns1;
  dBNode node;

  if(num_nodes == 0) return;

  node = forward ? nodes[0] : db_node_reverse(nodes[num_nodes-1]);
  kolist = kograph_get(kograph, node.key);
  num_occur = kograph_num(kograph, node.key);

  if(korun->len + num_occur == 0) { kmer_run_buf_reset(korun); return; }

  size_t max_paths = korun->len + num_occur;

  // Split korun buffer into two fixed size arrays
  // Extend current runs, and start new ones
  kmer_run_buf_capacity(korun, 2*max_paths);
  runs0 = korun->data;
  runs1 = korun->data + max_paths;
  nruns0 = korun->len;
  nruns1 = korun_extend(runs0, nruns0, node, kolist, num_occur, runs1,
                        pickup_at_first_node, qoffset);

  // Store runs that could not be extended
  for(i = 0; i < nruns0; i++) {
    if(!runs0[i].used && korun_len(runs0[i]) > min_len) {
      runs_ended->data[runs_ended->len++] = runs0[i];
    }
  }

  SWAP( runs0,  runs1);
  SWAP(nruns0, nruns1);

  for(i = 1; i < num_nodes && nruns0 > 0; i++)
  {
    node = forward ? nodes[i] : db_node_reverse(nodes[num_nodes-i-1]);
    kolist = kograph_get(kograph, node.key);
    num_occur = kograph_num(kograph, node.key);
    nruns1 = korun_extend(runs0, nruns0, node, kolist, num_occur, runs1, false, 0);

    // Store runs than have ended
    if(nruns0 < nruns1) {
      for(j = 0; j < nruns0; j++) {
        if(!runs0[j].used && korun_len(runs0[j]) > min_len) {
          runs_ended->data[runs_ended->len++] = runs0[j];
        }
      }
    }

    SWAP( runs0,  runs1);
    SWAP(nruns0, nruns1);
  }

  if(runs0 != korun->data) {
    // Copy final array into start position
    memmove(korun->data, runs0, nruns0 * sizeof(KOccurRun));
  }

  korun->len = nruns0;
}

