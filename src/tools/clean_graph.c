#include "global.h"
#include "util.h"
#include "file_util.h"
#include "supernode.h"
#include "prune_nodes.h"
#include "clean_graph.h"

#define DUMP_COVG_ARRSIZE 1000
#define DUMP_LEN_ARRSIZE 1000

// Define a vector of Covg
#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(covg_buf,CovgBuffer,Covg);


typedef struct
{
  const size_t nthreads, covg_threshold, min_keep_tip;
  CovgBuffer *cbufs;
  // uint64_t *covg_hist;
  uint64_t *covg_hist_init, *covg_hist_cleaned;
  uint64_t *covg_kmers_hist_init, *covg_kmers_hist_cleaned;
  uint64_t *len_hist_init, *len_hist_cleaned;
  const size_t covg_arrlen, len_arrlen;
  uint8_t *keep_flags;
  uint64_t num_tips,      num_low_covg_snodes,      num_tip_and_low_snodes;
  uint64_t num_tip_kmers, num_low_covg_snode_kmers, num_tip_and_low_snode_kmers;
  const dBGraph *db_graph;
} SupernodeCleaner;

// #define supernode_covg(covgs,len) supernode_covg_mean(covgs,len)
#define supernode_covg(covgs,len) supernode_read_starts(covgs,len)

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
static inline void fetch_coverages(dBNodeBuffer nbuf, CovgBuffer *cbuf,
                                   const dBGraph *db_graph)
{
  size_t i;
  covg_buf_reset(cbuf);
  covg_buf_capacity(cbuf, nbuf.len);
  cbuf->len = nbuf.len;
  for(i = 0; i < nbuf.len; i++)
    cbuf->data[i] = db_graph->col_covgs[nbuf.data[i].key];
}

static inline bool nodes_are_tip(dBNodeBuffer nbuf, const dBGraph *db_graph)
{
  Edges first = db_node_get_edges_union(db_graph, nbuf.data[0].key);
  Edges last = db_node_get_edges_union(db_graph, nbuf.data[nbuf.len-1].key);
  int in = edges_get_indegree(first, nbuf.data[0].orient);
  int out = edges_get_outdegree(last, nbuf.data[nbuf.len-1].orient);
  return (in+out <= 1);
}

static inline bool nodes_are_removable_tip(dBNodeBuffer nbuf,
                                           size_t min_keep_tip,
                                           const dBGraph *db_graph)
{
  return (nbuf.len < min_keep_tip && nodes_are_tip(nbuf, db_graph));
}


static void supernode_cleaner_alloc(SupernodeCleaner *cl, size_t nthreads,
                                    size_t covg_threshold, size_t min_keep_tip,
                                    uint8_t *keep_flags,
                                    const dBGraph *db_graph)
{
  size_t i;
  CovgBuffer *cbufs = ctx_calloc(nthreads, sizeof(CovgBuffer));
  for(i = 0; i < nthreads; i++)
    covg_buf_alloc(&cbufs[i], 1024);

  uint64_t *covg_hist_init, *covg_hist_cleaned;
  uint64_t *covg_kmers_hist_init, *covg_kmers_hist_cleaned;
  uint64_t *len_hist_init, *len_hist_cleaned;

  covg_hist_init          = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  covg_hist_cleaned       = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  covg_kmers_hist_init    = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  covg_kmers_hist_cleaned = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  len_hist_init           = ctx_calloc(DUMP_LEN_ARRSIZE, sizeof(uint64_t));
  len_hist_cleaned        = ctx_calloc(DUMP_LEN_ARRSIZE, sizeof(uint64_t));

  SupernodeCleaner tmp = {.nthreads = nthreads,
                          .covg_threshold = covg_threshold,
                          .min_keep_tip = min_keep_tip,
                          .cbufs = cbufs,
                          .covg_hist_init    = covg_hist_init,
                          .covg_hist_cleaned = covg_hist_cleaned,
                          .covg_kmers_hist_init = covg_kmers_hist_init,
                          .covg_kmers_hist_cleaned = covg_kmers_hist_cleaned,
                          .len_hist_init     = len_hist_init,
                          .len_hist_cleaned  = len_hist_cleaned,
                          .covg_arrlen = DUMP_COVG_ARRSIZE,
                          .len_arrlen = DUMP_LEN_ARRSIZE,
                          .keep_flags = keep_flags,
                          .num_tips = 0,
                          .num_low_covg_snodes = 0,
                          .num_tip_and_low_snodes = 0,
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
  ctx_free(cl->covg_hist_init);
  ctx_free(cl->covg_hist_cleaned);
  ctx_free(cl->covg_kmers_hist_init);
  ctx_free(cl->covg_kmers_hist_cleaned);
  ctx_free(cl->len_hist_init);
  ctx_free(cl->len_hist_cleaned);
  memset(cl, 0, sizeof(SupernodeCleaner));
}

static inline void supernode_get_covg(dBNodeBuffer nbuf, size_t threadid,
                                      void *arg)
{
  const SupernodeCleaner *cl = (const SupernodeCleaner*)arg;
  size_t covg, len;

  // Get coverage
  CovgBuffer *cbuf = &cl->cbufs[threadid];
  fetch_coverages(nbuf, cbuf, cl->db_graph);
  covg = supernode_covg(cbuf->data, cbuf->len);
  covg = MIN2(covg, cl->covg_arrlen-1);
  len = MIN2(cbuf->len, cl->len_arrlen-1);

  // Add to histograms
  __sync_fetch_and_add((volatile uint64_t *)&cl->covg_hist_init[covg], 1);
  __sync_fetch_and_add((volatile uint64_t *)&cl->covg_kmers_hist_init[covg], cbuf->len);
  __sync_fetch_and_add((volatile uint64_t *)&cl->len_hist_init[len], 1);
}

// Get coverage threshold for removing supernodes
// If `min_keep_tip` is > 0, tips shorter than `min_keep_tip` are not used
// in measuring supernode coverage.
// `visited`, should each be at least db_graph.ht.capcity bits long
//   and initialised to zero. On return, it will be 1 at each original kmer index
Covg cleaning_get_threshold(size_t num_threads, double seq_depth,
                            const char *covgs_csv_path,
                            const char *lens_csv_path,
                            uint8_t *visited, dBGraph *db_graph)
{
  // Estimate optimum cleaning threshold
  status("[cleaning] Calculating supernode statistiscs with %zu threads...",
         num_threads);

  // Get supernode coverages and lengths
  SupernodeCleaner cl;
  supernode_cleaner_alloc(&cl, num_threads, 0, 0, NULL, db_graph);

  supernodes_iterate(num_threads, visited, db_graph, supernode_get_covg, &cl);

  // Wipe visited kmer memory
  memset(visited, 0, roundup_bits2bytes(db_graph->ht.capacity));

  if(covgs_csv_path != NULL) {
    cleaning_write_covg_histogram(covgs_csv_path, cl.covg_hist_init,
                                  cl.covg_kmers_hist_init, cl.covg_arrlen);
  }

  if(lens_csv_path != NULL) {
    cleaning_write_len_histogram(lens_csv_path, cl.len_hist_init, cl.len_arrlen,
                                 db_graph->kmer_size);
  }

  // set threshold using histogram and genome size
  size_t threshold_est = cleaning_supernode_threshold(cl.covg_hist_init,
                                                      cl.covg_arrlen,
                                                      seq_depth,
                                                      db_graph);

  status("[cleaning] Recommended supernode cleaning threshold: < %zu",
         threshold_est);

  supernode_cleaner_dealloc(&cl);

  return threshold_est;
}

static inline void supernode_mark(dBNodeBuffer nbuf, size_t threadid,
                                  void *arg)
{
  SupernodeCleaner *cl = (SupernodeCleaner*)arg;
  bool low_covg_snode = false, removable_tip = false;
  size_t i, covg, len;

  CovgBuffer *cbuf = &cl->cbufs[threadid];
  fetch_coverages(nbuf, cbuf, cl->db_graph);
  covg = supernode_covg(cbuf->data, cbuf->len);
  low_covg_snode = (covg < cl->covg_threshold);

  // Remove tips
  removable_tip = nodes_are_removable_tip(nbuf, cl->min_keep_tip, cl->db_graph);

  if(low_covg_snode && removable_tip) {
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_tip_and_low_snodes, 1);
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_tip_and_low_snode_kmers, nbuf.len);
  } else if(low_covg_snode) {
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_low_covg_snodes, 1);
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_low_covg_snode_kmers, nbuf.len);
  } else if(removable_tip) {
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_tips, 1);
    __sync_fetch_and_add((volatile uint64_t *)&cl->num_tip_kmers, nbuf.len);
  } else {
    for(i = 0; i < nbuf.len; i ++)
      (void)bitset_set_mt(cl->keep_flags, nbuf.data[i].key);

    // Add to histograms
    covg = MIN2(covg, cl->covg_arrlen-1);
    len = MIN2(nbuf.len, cl->covg_arrlen-1);

    __sync_fetch_and_add((volatile uint64_t *)&cl->covg_hist_cleaned[covg], 1);
    __sync_fetch_and_add((volatile uint64_t *)&cl->covg_kmers_hist_cleaned[covg], cbuf->len);
    __sync_fetch_and_add((volatile uint64_t *)&cl->len_hist_cleaned[len], 1);
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
                 const char *covgs_csv_path, const char *lens_csv_path,
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

  if(covg_threshold > 0)
    status("[cleaning] Removing supernodes with coverage < %zu...", covg_threshold);
  if(min_keep_tip > 0)
    status("[cleaning] Removing tips shorter than %zu...", min_keep_tip);

  status("[cleaning]   using %zu threads", num_threads);

  // Mark nodes to keep
  SupernodeCleaner cl;
  supernode_cleaner_alloc(&cl, num_threads, covg_threshold, min_keep_tip,
                          keep, db_graph);
  supernodes_iterate(num_threads, visited, db_graph, supernode_mark, &cl);

  // Print numbers of kmers that are being removed

  char num_snodes_str[50], num_tips_str[50], num_tip_snodes_str[50];
  char num_snode_kmers_str[50], num_tip_kmers_str[50], num_tip_snode_kmers_str[50];
  ulong_to_str(cl.num_low_covg_snodes, num_snodes_str);
  ulong_to_str(cl.num_tips, num_tips_str);
  ulong_to_str(cl.num_tip_and_low_snodes, num_tip_snodes_str);
  ulong_to_str(cl.num_low_covg_snode_kmers, num_snode_kmers_str);
  ulong_to_str(cl.num_tip_kmers, num_tip_kmers_str);
  ulong_to_str(cl.num_tip_and_low_snode_kmers, num_tip_snode_kmers_str);

  status("[cleaning] Removing %s low coverage supernode%s [%s kmer%s], "
         "%s supernode tip%s [%s kmer%s] "
         "and %s of both [%s kmer%s]",
         num_snodes_str, util_plural_str(cl.num_low_covg_snodes),
         num_snode_kmers_str, util_plural_str(cl.num_low_covg_snode_kmers),
         num_tips_str, util_plural_str(cl.num_tips),
         num_tip_kmers_str, util_plural_str(cl.num_tip_kmers),
         num_tip_snodes_str,
         num_tip_snode_kmers_str, util_plural_str(cl.num_tip_and_low_snode_kmers));

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

  if(covgs_csv_path != NULL) {
    cleaning_write_covg_histogram(covgs_csv_path, cl.covg_hist_cleaned,
                                  cl.covg_kmers_hist_cleaned, cl.covg_arrlen);
  }

  if(lens_csv_path != NULL) {
    cleaning_write_len_histogram(lens_csv_path, cl.len_hist_cleaned,
                                 cl.len_arrlen, db_graph->kmer_size);
  }

  supernode_cleaner_dealloc(&cl);
}

static FILE* _open_histogram_file(const char *path, const char *name)
{
  FILE *fout;
  status("[cleaning] Writing %s distribution to: %s", name, futil_outpath_str(path));

  if(strcmp(path,"-") == 0) return stdout;
  if((fout = fopen(path, "w")) == NULL)
    warn("Couldn't write %s distribution to file: %s", name, path);
  return fout;
}

void cleaning_write_covg_histogram(const char *path,
                                   const uint64_t *covg_hist,
                                   const uint64_t *kmer_hist,
                                   size_t len)
{
  ctx_assert(len >= 2);
  ctx_assert(covg_hist[0] == 0);
  ctx_assert(kmer_hist[0] == 0);
  size_t i, end;

  FILE *fout = _open_histogram_file(path, "supernode coverage");
  if(fout == NULL) return;

  fprintf(fout, "Covg,NumSupernodes,NumKmers\n");
  for(end = len-1; end > 1 && covg_hist[end] == 0; end--) {}
  fprintf(fout, "1,%"PRIu64",%"PRIu64"\n", covg_hist[1], kmer_hist[1]);
  for(i = 2; i <= end; i++) {
    if(covg_hist[i] > 0)
      fprintf(fout, "%zu,%"PRIu64",%"PRIu64"\n", i, covg_hist[i], kmer_hist[i]);
  }
  fclose(fout);
}

void cleaning_write_len_histogram(const char *path,
                                  const uint64_t *hist, size_t len,
                                  size_t kmer_size)
{
  ctx_assert(len >= 2);
  ctx_assert(hist[0] == 0);
  size_t i, end;

  FILE *fout = _open_histogram_file(path, "supernode length");
  if(fout == NULL) return;

  fprintf(fout, "SupernodeKmerLength,bp,Count\n");
  for(end = len-1; end > 1 && hist[end] == 0; end--) {}
  fprintf(fout, "1,%zu,%"PRIu64"\n", kmer_size, hist[1]);
  for(i = 2; i <= end; i++) {
    if(hist[i] > 0)
      fprintf(fout, "%zu,%zu,%"PRIu64"\n", i, kmer_size+i-1, hist[i]);
  }
  fclose(fout);
}
