#include "global.h"
#include "util.h"
#include "supernode.h"
#include "prune_nodes.h"
#include "clean_graph.h"

#define DUMP_COVG_ARRSIZE 1000

// Define a vector of Covg
#include "objbuf_macro.h"
create_objbuf(covg_buf,CovgBuffer,Covg);

#define supernode_covg(covgs,len) supernode_covg_mean(covgs,len)
// #define supernode_covg(covgs,len) supernode_read_starts(covgs,len)

static inline size_t supernode_covg_mean(Covg *covgs, size_t len)
{
  ctx_assert(len > 0);
  size_t i, sum = 0;
  for(i = 0; i < len; i++) sum += covgs[i];
  return (sum+len/2) / len; // round to nearest integer
}

// Calculate cleaning threshold for supernodes from a given distribution
// of supernode coverages
size_t cleaning_supernode_threshold(const uint64_t *covgs, size_t len,
                                    double seq_depth,
                                    const dBGraph *db_graph)
{
  ctx_assert(len > 5);
  ctx_assert(db_graph->ht.num_kmers > 0);

  size_t i, d1len = len-2, d2len = len-3, f1, f2;
  double *tmp = ctx_malloc((d1len+d2len) * sizeof(double));
  double *delta1 = tmp, *delta2 = tmp + d1len;

  // Get sequencing depth from coverage
  uint64_t covg_sum = 0, capacity = db_graph->ht.capacity * db_graph->num_of_cols;
  for(i = 0; i < capacity; i++) covg_sum += db_graph->col_covgs[i];
  double seq_depth_est = (double)covg_sum / db_graph->ht.num_kmers;

  status("[cleaning] Kmer depth before cleaning supernodes: %.2f", seq_depth_est);
  if(seq_depth <= 0) seq_depth = seq_depth_est;
  else status("[cleaning] Using sequence depth argument: %f", seq_depth);

  size_t fallback_thresh = (size_t)MAX2(1, (seq_depth+1)/2);

  // +1 to ensure covgs is never 0
  for(i = 0; i < d1len; i++) delta1[i] = (double)(covgs[i+1]+1) / (covgs[i+2]+1);

  d1len = i;
  d2len = d1len - 1;

  if(d1len <= 2) {
    status("[cleaning]  (using fallback1)\n");
    ctx_free(tmp);
    return fallback_thresh;
  }

  // d2len is d1len-1
  for(i = 0; i < d2len; i++) delta2[i] = delta1[i] / delta1[i+1];

  for(f1 = 0; f1 < d1len && delta1[f1] >= 1; f1++);
  for(f2 = 0; f2 < d2len && delta2[f2] > 1; f2++);

  ctx_free(tmp);

  if(f1 < d1len && f1 < (seq_depth*0.75))
  { status("[cleaning]   (using f1)"); return f1+1; }
  else if(f2 < d2len)
  { status("[cleaning]   (using f2)"); return f2+1; }
  else
  { status("[cleaning]   (using fallback1)"); return fallback_thresh+1; }
}

// Get coverages from nodes in nbuf, store in cbuf
static inline void fetch_coverages(const dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                                   const dBGraph *db_graph)
{
  size_t i;
  covg_buf_reset(cbuf);
  covg_buf_ensure_capacity(cbuf, nbuf->len);
  cbuf->len = nbuf->len;
  for(i = 0; i < nbuf->len; i++)
    cbuf->data[i] = db_graph->col_covgs[nbuf->data[i].key];
}

static inline bool nodes_are_tip(const dBNodeBuffer *nbuf,
                                 const dBGraph *db_graph)
{
  Edges first = db_node_get_edges_union(db_graph, nbuf->data[0].key);
  Edges last = db_node_get_edges_union(db_graph, nbuf->data[nbuf->len-1].key);
  int in = edges_get_indegree(first, nbuf->data[0].orient);
  int out = edges_get_outdegree(last, nbuf->data[nbuf->len-1].orient);
  return (in+out <= 1);
}

static inline bool nodes_are_removable_tip(const dBNodeBuffer *nbuf,
                                           size_t min_keep_tip,
                                           const dBGraph *db_graph)
{
  return (nbuf->len < min_keep_tip && nodes_are_tip(nbuf, db_graph));
}



typedef struct
{
  const size_t nthreads, covg_threshold, min_keep_tip;
  CovgBuffer *cbufs;
  uint64_t *covg_hist;
  const size_t covg_arrlen;
  uint8_t *keep_flags;
  uint64_t num_tip_kmers, num_low_covg_snode_kmers, num_tip_and_low_snode_kmers;
  const dBGraph *db_graph;
} SupernodeCleaner;

static void supernode_cleaner_alloc(SupernodeCleaner *cl, size_t nthreads,
                                    size_t covg_threshold, size_t min_keep_tip,
                                    uint64_t *covg_hist, size_t covg_arrlen,
                                    uint8_t *keep_flags,
                                    const dBGraph *db_graph)
{
  size_t i;
  CovgBuffer *cbufs = ctx_calloc(nthreads, sizeof(CovgBuffer));
  for(i = 0; i < nthreads; i++)
    covg_buf_alloc(&cbufs[i], 1024);

  SupernodeCleaner tmp = {.nthreads = nthreads,
                          .covg_threshold = covg_threshold,
                          .min_keep_tip = min_keep_tip,
                          .cbufs = cbufs,
                          .covg_hist = covg_hist,
                          .covg_arrlen = covg_arrlen,
                          .keep_flags = keep_flags,
                          .num_tip_kmers = 0,
                          .num_low_covg_snode_kmers = 0,
                          .num_tip_and_low_snode_kmers = 0,
                          .db_graph = db_graph};

  memcpy(cl, &tmp, sizeof(SupernodeCleaner));
}

static void supernode_cleaner_dealloc(SupernodeCleaner *cl)
{
  size_t i;
  for(i = 0; i < cl->nthreads; i++)
    covg_buf_dealloc(&cl->cbufs[i]);
  ctx_free(cl->cbufs);
  memset(cl, 0, sizeof(SupernodeCleaner));
}

static inline void supernode_get_covg(const dBNodeBuffer *nbuf, size_t threadid,
                                      void *arg)
{
  const SupernodeCleaner *cl = (const SupernodeCleaner*)arg;

  // Check if this is not a tip that we should ignore
  if(!nodes_are_removable_tip(nbuf, cl->min_keep_tip, cl->db_graph))
  {
    // Get coverage
    CovgBuffer *cbuf = &cl->cbufs[threadid];
    fetch_coverages(nbuf, cbuf, cl->db_graph);
    size_t reads_arriving = supernode_covg(cbuf->data, cbuf->len);
    reads_arriving = MIN2(reads_arriving, cl->covg_arrlen-1);

    // Add to histogram
    __sync_fetch_and_add((volatile uint64_t *)&cl->covg_hist[reads_arriving], 1);
  }
}

// Get coverage threshold for removing supernodes
// If `min_keep_tip` is > 0, tips shorter than `min_keep_tip` are not used
// in measuring supernode coverage.
// `visited`, should each be at least db_graph.ht.capcity bits long
//   and initialised to zero. On return, it will be 1 at each original kmer index
Covg cleaning_get_threshold(size_t num_threads, size_t min_keep_tip,
                            double seq_depth, const char *dump_covgs,
                            uint8_t *visited, dBGraph *db_graph)
{
  // Estimate optimum cleaning threshold
  status("[cleaning] Calculating supernode cleaning threshold with %zu threads...",
         num_threads);

  // Get supernode coverages
  uint64_t *covg_hist = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));

  SupernodeCleaner cleaner;
  supernode_cleaner_alloc(&cleaner, num_threads, 0, min_keep_tip,
                          covg_hist, DUMP_COVG_ARRSIZE, NULL, db_graph);

  supernodes_iterate(num_threads, visited, db_graph,
                     supernode_get_covg, &cleaner);

  supernode_cleaner_dealloc(&cleaner);

  // Wipe memory
  memset(visited, 0, roundup_bits2bytes(db_graph->ht.capacity));

  if(dump_covgs != NULL)
    cleaning_dump_covg_histogram(dump_covgs, covg_hist, DUMP_COVG_ARRSIZE);

  // set threshold using histogram and genome size
  size_t threshold_est = cleaning_supernode_threshold(covg_hist,
                                                      DUMP_COVG_ARRSIZE,
                                                      seq_depth,
                                                      db_graph);

  status("[cleaning] Recommended supernode cleaning threshold: < %zu",
         threshold_est);

  ctx_free(covg_hist);

  return threshold_est;
}

static inline void supernode_mark(const dBNodeBuffer *nbuf, size_t threadid,
                                   void *arg)
{
  const SupernodeCleaner *cl = (const SupernodeCleaner*)arg;
  bool low_covg_snode = false, removable_tip = false;
  size_t i;

  if(cl->covg_threshold > 0)
  {
    // Clean supernodes by coverage
    CovgBuffer *cbuf = &cl->cbufs[threadid];
    fetch_coverages(nbuf, cbuf, cl->db_graph);
    size_t reads_arriving = supernode_covg(cbuf->data, cbuf->len);
    low_covg_snode = (reads_arriving < cl->covg_threshold);
  }

  // Remove tips
  if(nodes_are_removable_tip(nbuf, cl->min_keep_tip, cl->db_graph))
  {
    removable_tip = true;
  }

  if(low_covg_snode && removable_tip)
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_tip_and_low_snode_kmers, nbuf->len);
  else if(low_covg_snode)
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_low_covg_snode_kmers, nbuf->len);
  else if(removable_tip)
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_tip_kmers, nbuf->len);
  else {
    for(i = 0; i < nbuf->len; i ++)
      bitset_set_mt(cl->keep_flags, nbuf->data[i].key);
  }
}

// Remove low coverage supernodes and clip tips
// - Remove supernodes with coverage < `covg_threshold`
// - Remove tips shorter than `min_keep_tip`
// `visited`, `keep` should each be at least db_graph.ht.capcity bits long
//   and initialised to zero. On return,
//   `visited` will be 1 at each original kmer index
//   `keep` will be 1 at each retained kmer index
void clean_graph(size_t num_threads, size_t covg_threshold, size_t min_keep_tip,
                 uint8_t *visited, uint8_t *keep, dBGraph *db_graph)
{
  ctx_assert(db_graph->num_of_cols == 1);
  ctx_assert(db_graph->num_edge_cols > 0);

  size_t init_nkmers = db_graph->ht.num_kmers;

  if(db_graph->ht.num_kmers == 0) return;
  if(covg_threshold == 0 && min_keep_tip == 0) {
    warn("[cleaning] No cleaning specified");
    return;
  }

  if(covg_threshold > 0 && min_keep_tip > 0)
    status("[cleaning] Removing supernodes with coverage < %zu and "
           "tips shorter than %zu...", covg_threshold, min_keep_tip);
  else if(covg_threshold > 0)
    status("[cleaning] Removing supernodes with coverage < %zu", covg_threshold);
  else
    status("[cleaning] Removing tips shorter than %zu...", min_keep_tip);

  status("[cleaning]   using %zu threads", num_threads);

  // Mark nodes to keep
  SupernodeCleaner cleaner;
  supernode_cleaner_alloc(&cleaner, num_threads, covg_threshold, min_keep_tip,
                          NULL, 0, keep, db_graph);
  supernodes_iterate(num_threads, visited, db_graph, supernode_mark, &cleaner);

  status("[cleaning] Removing %zu supernode kmers, %zu tip kmers and %zu of both",
         (size_t)cleaner.num_low_covg_snode_kmers, (size_t)cleaner.num_tip_kmers,
         (size_t)cleaner.num_tip_and_low_snode_kmers);

  supernode_cleaner_dealloc(&cleaner);

  // Remove nodes not marked to keep
  prune_nodes_lacking_flag(num_threads, keep, db_graph);

  // Wipe memory
  memset(visited, 0, roundup_bits2bytes(db_graph->ht.capacity));
  memset(keep, 0, roundup_bits2bytes(db_graph->ht.capacity));

  // Print status update
  char remain_nkmers_str[100], removed_nkmers_str[100];
  size_t remain_nkmers = db_graph->ht.num_kmers;
  size_t removed_nkmers = init_nkmers - remain_nkmers;
  ulong_to_str(remain_nkmers, remain_nkmers_str);
  ulong_to_str(removed_nkmers, removed_nkmers_str);
  status("[cleaning] Remaining kmers: %s removed: %s (%.1f%%)",
         remain_nkmers_str, removed_nkmers_str,
         (100.0*removed_nkmers)/init_nkmers);
}

void cleaning_dump_covg_histogram(const char *path, uint64_t *hist, size_t len)
{
  if(len == 0) return;
  ctx_assert(len >= 2);
  ctx_assert(hist[0] == 0);

  status("[cleaning] Writing covg distribution to: %s", path);
  size_t i, end;
  FILE *fout = fopen(path, "w");

  if(fout == NULL) {
    warn("Couldn't write size distribution! file: %s", path);
    return;
  }

  fprintf(fout, "Covg,Supernodes\n");
  fprintf(fout, "1,%zu\n", (size_t)hist[1]);
  for(end = len-1; end > 1 && hist[end] == 0; end--);
  for(i = 2; i <= end; i++) {
    if(hist[i] > 0)
      fprintf(fout, "%zu,%zu\n", i, (size_t)hist[i]);
  }
  fclose(fout);
}
