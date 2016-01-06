#include "global.h"
#include "util.h"
#include "file_util.h"
#include "supernode.h"
#include "prune_nodes.h"
#include "clean_graph.h"

#include "carrays/carrays.h" // gca_median()

#include <math.h> // lgamma, tgamma
#include <float.h> // DBL_MAX

#define DUMP_COVG_ARRSIZE 1000
#define DUMP_LEN_ARRSIZE 1000

// Define a vector of Covg
#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(covg_buf,CovgBuffer,Covg);

// Find cutoff by finding first coverage level where errors make up less than
// `fdr` of total coverage
// returns -1 if not found
static inline int pick_cutoff_with_fdr_thresh(const double *e_covg,
                                              const uint64_t *kmer_covg,
                                              size_t arrlen, double fdr)
{
  size_t i;
  for(i = 1; i < arrlen; i++) {
    // printf(" %zu: %f %zu test: %f < %f\n", i, e_covg[i], kmer_covg[i],
    //                                       e_covg[i] / kmer_covg[i], fdr);
    if(e_covg[i] / kmer_covg[i] <= fdr) {
      return i;
    }
  }
  return -1;
}

// Get highest cutoff where false-positives < false-negatives
// i.e. proportion of real kmers we are removing is less than the
//      proportion of bad kmers we are keeping
// returns -1 if not found
static inline int pick_cutoff_FP_lt_FN(const double *e_covg, double e_total,
                                      const uint64_t *kmer_covg, uint64_t d_total,
                                      size_t arrlen)
{
  size_t i;
  // for(i = 1; i < arrlen; i++) { printf("  %zu", kmer_covg[i]); } printf("\n");
  // for(i = 1; i < arrlen; i++) { printf("  %f", e_covg[i]); } printf("\n");
  // printf(" e_total: %f d_total: %f\n", e_total, (double)d_total);
  double e_rem = e_total, d_rem = d_total;
  double e_sum = 0, d_sum = 0;
  for(i = 1; i < arrlen; i++) {
    e_sum += e_covg[i];
    d_sum += kmer_covg[i];
    e_rem -= e_covg[i];
    d_rem -= kmer_covg[i];
    // printf(" %zu: e_total: %f d_total: %f\n", i, e_rem, d_rem);
    if(1-e_sum/d_sum > e_rem/d_rem) {
      return i;
    }
  }
  return -1;
}

static inline int pick_cutoff_loss_vs_error(const double *e_covg,
                                            double e_total,
                                            const uint64_t *kmer_covg,
                                            size_t arrlen)
{
  size_t i;
  // for(i = 1; i < arrlen; i++) { printf("  %zu", kmer_covg[i]); } printf("\n");
  // for(i = 1; i < arrlen; i++) { printf("  %f", e_covg[i]); } printf("\n");
  // printf(" e_total: %f d_total: %f\n", e_total, (double)d_total);
  double e_rem = e_total;
  double e_sum = 0, d_sum = 0;
  for(i = 1; i < arrlen; i++) {
    e_sum += e_covg[i];
    d_sum += kmer_covg[i];
    e_rem -= e_covg[i];
    double lost_seq = (d_sum-e_sum);
    double rem_err = e_rem;
    // printf(" %zu: e_total: %f d_total: %f\n", i, e_rem, d_rem);
    if(lost_seq > rem_err) return i;
  }
  return -1;
}

static inline void cutoff_get_FP_FN(const double *e_covg, double e_total,
                                    const uint64_t *kmer_covg, uint64_t d_total,
                                    size_t cutoff,
                                    double *false_pos, double *false_neg)
{
  size_t i;
  double e_rem = e_total, d_rem = d_total;
  double e_sum = 0, d_sum = 0;
  for(i = 1; i < cutoff; i++) {
    e_sum += e_covg[i];
    d_sum += kmer_covg[i];
    e_rem -= e_covg[i];
    d_rem -= kmer_covg[i];
  }
  *false_pos = 1-e_sum/d_sum;
  *false_neg = e_rem/d_rem;
}

// Check if at least `frac_covg_kept` coverage is kept when using threshold
static inline bool is_cutoff_good(const uint64_t *kmer_covg, size_t arrlen,
                                  size_t cutoff, double frac_covg_kept)
{
  uint64_t kmers_below = 0, kmers_above = 0;
  size_t i;
  for(i = 0;      i < cutoff; i++) kmers_below += kmer_covg[i]*i;
  for(i = cutoff; i < arrlen; i++) kmers_above += kmer_covg[i]*i;

  // At least 20% of kmers should be kept
  return !arrlen || // any cutoff is good if no kmers
         ((double)kmers_above/(kmers_below+kmers_above) >= frac_covg_kept);
}

/**
 * Pick a cleaning threshold from kmer coverage histogram. Assumes low coverage
 * kmers are all due to error. Fits a poisson with a gamma distributed mean.
 * Then chooses a cleaning threshold such than FDR (uncleaned kmers) occur at a
 * rate of < the FDR paramater.
 *
 * Translated from Gil McVean's initial proposed method in R code
 *
 * @param kmer_covg Histogram of kmer counts at coverages 1,2,.. arrlen-1
 * @param arrlen    Length of array kmer_covg
 * @param alpha_est_ptr If not NULL, used to return estimate for alpha
 * @param beta_est_ptr  If not NULL, used to return estimate for beta
 * @return -1 if no cut-off satisfies FDR, otherwise returns coverage cutoff
 */
int cleaning_pick_kmer_threshold(const uint64_t *kmer_covg, size_t arrlen,
                                 double *alpha_est_ptr, double *beta_est_ptr,
                                 double *false_pos_ptr, double *false_neg_ptr)
{
  ctx_assert(arrlen >= 10);
  ctx_assert2(kmer_covg[0] == 0, "Shouldn't see any kmers with coverage zero");

  size_t i, min_a_est_idx = 0;
  double r1, r2, rr, min_a_est = DBL_MAX, tmp;
  double aa, faa, a_est, b_est, c0;

  r1 = (double)kmer_covg[2] / kmer_covg[1];
  r2 = (double)kmer_covg[3] / kmer_covg[2];
  rr = r2 / r1;

  // printf("r1: %.2f r2: %.2f rr: %.2f\n", r1, r2, rr);

  // iterate aa = { 0.01, 0.02, ..., 1.99, 2.00 }
  // find aa value that minimises abs(faa-rr)
  for(i = 1; i <= 200; i++)
  {
    aa = i*0.01;
    faa = tgamma(aa)*tgamma(aa+2) / (2*pow(tgamma(aa+1),2));
    tmp = fabs(faa-rr);
    if(tmp < min_a_est) { min_a_est = tmp; min_a_est_idx = i; }
  }

  // a_est, b_est are estimates for alpha, beta of gamma distribution
  a_est = min_a_est_idx*0.01;
  b_est = tgamma(a_est + 1.0) / (r1 * tgamma(a_est)) - 1.0;
  b_est = MAX2(b_est, 1); // Avoid beta values <1
  c0 = kmer_covg[1] * pow(b_est/(1+b_est),-a_est);

  if(alpha_est_ptr) *alpha_est_ptr = a_est;
  if(beta_est_ptr)  *beta_est_ptr  = b_est;

  // printf("min_a_est_idx: %zu\n", min_a_est_idx);
  // printf("a_est: %f b_est %f c0: %f\n", a_est, b_est, c0);

  // keep coverage estimates on the stack - this should be ok
  double e_covg_tmp, e_covg[arrlen];
  double e_total = 0;
  uint64_t d_total = 0;

  // Calculate some values here for speed
  double log_b_est          = log(b_est);
  double log_one_plus_b_est = log(1 + b_est);
  double lgamma_a_est       = lgamma(a_est);

  // note: lfactorial(x) = lgamma(x+1)

  for(i = 1; i < arrlen; i++)
  {
    e_covg_tmp = a_est * log_b_est - lgamma_a_est - lgamma(i)
                   + lgamma(a_est + i - 1)
                   - (a_est + i - 1) * log_one_plus_b_est;
    e_covg[i] = exp(e_covg_tmp) * c0;
    e_total += e_covg[i];
    d_total += kmer_covg[i];
  }

  // for(i = 1; i < MIN2(arrlen,100); i++)
  //   printf("  %zu: %f %zu\n", i, e_covg[i], (size_t)kmer_covg[i]);

  int cutoff = -1;

  // Find cutoff by finding first coverage level where errors make up less than
  // 0.1% of total coverage
  cutoff = pick_cutoff_with_fdr_thresh(e_covg, kmer_covg, arrlen, 0.001);
  // printf("A cutoff: %i\n", cutoff);

  // Pick highest cutoff that keeps FP < FN
  if(cutoff < 0)
    cutoff = pick_cutoff_FP_lt_FN(e_covg, e_total, kmer_covg, d_total, arrlen);

  if(cutoff < 0)
    cutoff = pick_cutoff_loss_vs_error(e_covg, e_total, kmer_covg, arrlen);

  // printf("B cutoff: %i\n", cutoff);

  if(cutoff < 0) return -1;

  // printf("C cutoff: %i\n", cutoff);

  // Check cutoff keeps at least 20% of coverage
  // (WGS should be much higher, Exome sequencing needs low cutoff)
  if(!is_cutoff_good(kmer_covg, arrlen, cutoff, 0.2)) return -1;

  // printf("D cutoff: %i\n", cutoff);

  // Calculate FP,FN rates
  if(false_pos_ptr || false_neg_ptr) {
    double false_pos = 0, false_neg = 0;
    cutoff_get_FP_FN(e_covg, e_total, kmer_covg, d_total, cutoff,
                     &false_pos, &false_neg);
    // printf("  FP: %f, FN: %f\n", false_pos, false_neg);
    if(false_pos_ptr) *false_pos_ptr = false_pos;
    if(false_neg_ptr) *false_neg_ptr = false_neg;
  }

  // printf(" kmers_above : %zu / (%zu + %zu) = %f\n",
  //        kmers_above, kmers_below, kmers_above,
  //        (double)kmers_above/(kmers_below+kmers_above));

  // printf("cutoff: %i\n", cutoff);

  // printf(" cutoff: %zu fdr: %f fdr_limit: %f good: %i\n",
  //        cutoff, fdr, fdr_limit, (int)good_cutoff);

  return cutoff;
}

typedef struct
{
  const size_t nthreads, covg_threshold, min_keep_tip;
  CovgBuffer *cbufs;
  uint64_t *kmer_covgs_init, *kmer_covgs_clean;
  uint64_t *unitig_covgs_init, *unitig_covg_clean;
  uint64_t *len_hist_init, *len_hist_clean;
  const size_t covg_arrsize, len_arrsize;
  uint8_t *keep_flags;
  uint64_t num_tips,      num_low_covg_snodes,      num_tip_and_low_snodes;
  uint64_t num_tip_kmers, num_low_covg_snode_kmers, num_tip_and_low_snode_kmers;
  const dBGraph *db_graph;
} UnitigCleaner;

// Get coverages from nodes in nbuf, store in cbuf
static inline void fetch_coverages(dBNodeBuffer nbuf, CovgBuffer *cbuf,
                                   const dBGraph *db_graph)
{
  ctx_assert(db_graph->num_of_cols == 1);
  size_t i;
  covg_buf_reset(cbuf);
  covg_buf_capacity(cbuf, nbuf.len);
  cbuf->len = nbuf.len;
  for(i = 0; i < nbuf.len; i++)
    cbuf->b[i] = db_graph->col_covgs[nbuf.b[i].key];
}

static inline bool nodes_are_tip(dBNodeBuffer nbuf, const dBGraph *db_graph)
{
  Edges first = db_node_get_edges_union(db_graph, nbuf.b[0].key);
  Edges last = db_node_get_edges_union(db_graph, nbuf.b[nbuf.len-1].key);
  int in = edges_get_indegree(first, nbuf.b[0].orient);
  int out = edges_get_outdegree(last, nbuf.b[nbuf.len-1].orient);
  return (in+out <= 1);
}

static inline bool nodes_are_removable_tip(dBNodeBuffer nbuf,
                                           size_t min_keep_tip,
                                           const dBGraph *db_graph)
{
  return (nbuf.len < min_keep_tip && nodes_are_tip(nbuf, db_graph));
}


static void unitig_cleaner_alloc(UnitigCleaner *cl, size_t nthreads,
                                 size_t covg_threshold, size_t min_keep_tip,
                                 uint8_t *keep_flags,
                                 const dBGraph *db_graph)
{
  size_t i;
  CovgBuffer *cbufs = ctx_calloc(nthreads, sizeof(CovgBuffer));
  for(i = 0; i < nthreads; i++)
    covg_buf_alloc(&cbufs[i], 1024);

  uint64_t *kmer_covgs_init, *kmer_covgs_clean;
  uint64_t *unitig_covgs_init, *unitig_covg_clean;
  uint64_t *len_hist_init, *len_hist_clean;

  kmer_covgs_init      = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  kmer_covgs_clean    = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  unitig_covgs_init    = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  unitig_covg_clean  = ctx_calloc(DUMP_COVG_ARRSIZE, sizeof(uint64_t));
  len_hist_init        = ctx_calloc(DUMP_LEN_ARRSIZE,  sizeof(uint64_t));
  len_hist_clean     = ctx_calloc(DUMP_LEN_ARRSIZE,  sizeof(uint64_t));

  UnitigCleaner tmp = {.nthreads = nthreads,
                       .covg_threshold = covg_threshold,
                       .min_keep_tip = min_keep_tip,
                       .cbufs = cbufs,
                       .kmer_covgs_init   = kmer_covgs_init,
                       .kmer_covgs_clean  = kmer_covgs_clean,
                       .unitig_covgs_init = unitig_covgs_init,
                       .unitig_covg_clean = unitig_covg_clean,
                       .covg_arrsize      = DUMP_COVG_ARRSIZE,
                       .len_hist_init   = len_hist_init,
                       .len_hist_clean  = len_hist_clean,
                       .len_arrsize     = DUMP_LEN_ARRSIZE,
                       .keep_flags = keep_flags,
                       .num_tips = 0,
                       .num_low_covg_snodes = 0,
                       .num_tip_and_low_snodes = 0,
                       .num_tip_kmers = 0,
                       .num_low_covg_snode_kmers = 0,
                       .num_tip_and_low_snode_kmers = 0,
                       .db_graph = db_graph};

  memcpy(cl, &tmp, sizeof(UnitigCleaner));
}

static void unitig_cleaner_dealloc(UnitigCleaner *cl)
{
  size_t i;
  for(i = 0; i < cl->nthreads; i++)
    covg_buf_dealloc(&cl->cbufs[i]);
  ctx_free(cl->cbufs);
  ctx_free(cl->kmer_covgs_init);
  ctx_free(cl->kmer_covgs_clean);
  ctx_free(cl->unitig_covgs_init);
  ctx_free(cl->unitig_covg_clean);
  ctx_free(cl->len_hist_init);
  ctx_free(cl->len_hist_clean);
  memset(cl, 0, sizeof(UnitigCleaner));
}

// Returns unitig coverage
static inline uint64_t update_kmer_covg_hist(uint64_t *kcovg_hist, size_t covgsize,
                                             uint64_t *ucovg_hist, size_t ucovgsize,
                                             uint64_t *len_hist, size_t lensize,
                                             CovgBuffer *cbuf)
{
  uint64_t i, kcovg, unitig_covg, len;

  // Histogram is of each kmer coverage
  for(i = 0; i < cbuf->len; i++) {
    kcovg = MIN2(cbuf->b[i], covgsize-1);
    __sync_fetch_and_add((volatile uint64_t *)&kcovg_hist[kcovg], 1);
  }

  // Length histgogram
  len = MIN2(cbuf->len, lensize-1);
  __sync_fetch_and_add((volatile uint64_t *)&len_hist[len], 1);

  // Mean covg histogram
  // size_t sum_covg;
  // for(i = sum_covg = 0; i < cbuf->len; i++) sum_covg += cbuf->b[i];
  // unitig_covg = sum_covg / cbuf->len;

  // Median coverage
  unitig_covg = gca_median_uint32(cbuf->b, cbuf->len);
  unitig_covg = MIN2(unitig_covg, ucovgsize-1);
  __sync_fetch_and_add((volatile uint64_t *)&ucovg_hist[unitig_covg], 1);

  return unitig_covg;
}

static inline void unitig_get_covg(dBNodeBuffer nbuf, size_t threadid, void *arg)
{
  const UnitigCleaner *cl = (const UnitigCleaner*)arg;

  // Load coverage into buffer
  CovgBuffer *cbuf = &cl->cbufs[threadid];
  fetch_coverages(nbuf, cbuf, cl->db_graph);

  // Update before-cleaning histograms
  update_kmer_covg_hist(cl->kmer_covgs_init, cl->covg_arrsize,
                        cl->unitig_covgs_init, cl->covg_arrsize,
                        cl->len_hist_init, cl->len_arrsize,
                        cbuf);
}

/*
//
// We could just iterate over kmers instead of unitigs to get coverage hist
// Currently we iterate over unitigs instead which is slower.
//
typedef struct {
  size_t nthreads;
  UnitigCleaner *cl;
} KmerCleanerIterator;

static inline void kmer_get_covg_node(hkey_t hkey, void *arg)
{
  const UnitigCleaner *cl = (const UnitigCleaner*)arg;
  size_t covg = db_node_sum_covg(cl->db_graph, hkey);
  covg = MIN2(covg, cl->covg_arrsize-1);
  __sync_fetch_and_add((volatile uint64_t *)&cl->kmer_covgs_init[covg], 1);
}

static void kmer_get_covg(void *arg, size_t threadid)
{
  const KmerCleanerIterator *kcl = (const KmerCleanerIterator*)arg;
  HASH_ITERATE_PART(&kcl->db_graph->ht, threadid, kcl->nthreads,
                    kmer_get_covg_node, kcl->cl);
}
*/

/**
 * Get coverage threshold for removing unitigs
 *
 * @param visited should be at least db_graph.ht.capcity bits long and initialised
 *                to zero. On return, it will be 1 at each original kmer index
 * @param covgs_csv_path
 * @param lens_csv_path  paths to files to write CSV histogram of unitigs
                         coverages and lengths BEFORE ANY CLEANING.
 *                       If NULL these are ignored.
 * @return threshold to clean or -1 on error
 */
int cleaning_get_threshold(size_t num_threads,
                           const char *covgs_csv_path,
                           const char *lens_csv_path,
                           uint8_t *visited,
                           const dBGraph *db_graph)
{
  // Estimate optimum cleaning threshold
  status("[cleaning] Calculating unitig stats with %zu threads...", num_threads);
  status("[cleaning]   Using kmer gamma method");

  // Get kmer coverages and unitig lengths
  UnitigCleaner cl;
  unitig_cleaner_alloc(&cl, num_threads, 0, 0, NULL, db_graph);
  supernodes_iterate(num_threads, visited, db_graph, unitig_get_covg, &cl);

  // Get kmer coverage only (faster)
  // KmerCleanerIterator kcls[nthreads];
  // for(i = 0; i < nthreads; i++)
  //   kcls[i] = (KmerCleanerIterator){.threadid = i, .nthreads = nthreads, .cl = &cl};
  // util_run_threads(kcls, nthreads, sizeof(kcls[0]), nthreads, kmer_get_covg);

  // Wipe visited kmer memory
  memset(visited, 0, roundup_bits2bytes(db_graph->ht.capacity));

  if(covgs_csv_path != NULL) {
    cleaning_write_covg_histogram(covgs_csv_path,
                                  cl.kmer_covgs_init,
                                  cl.unitig_covgs_init,
                                  cl.covg_arrsize);
  }

  if(lens_csv_path != NULL) {
    cleaning_write_len_histogram(lens_csv_path,
                                 cl.len_hist_init,
                                 cl.len_arrsize,
                                 db_graph->kmer_size);
  }

  // set threshold using histogram and genome size
  double alpha = 0, beta = 0, false_pos = 0, false_neg = 0;
  int threshold_est = cleaning_pick_kmer_threshold(cl.kmer_covgs_init,
                                                   cl.covg_arrsize,
                                                   &alpha, &beta,
                                                   &false_pos, &false_neg);

  if(threshold_est < 0)
    warn("Cannot pick a cleaning threshold");
  else {
    status("[cleaning] alpha=%f, beta=%f FP=%f FN=%f",
           alpha, beta, false_pos, false_neg);
    status("[cleaning] Recommended unitig cleaning threshold: < %i",
           threshold_est);
  }

  unitig_cleaner_dealloc(&cl);

  return threshold_est;
}

/**
 * Mark a unitig to keep or delete. Update stats on decision.
 */
static inline void unitig_mark(dBNodeBuffer nbuf, size_t threadid, void *arg)
{
  UnitigCleaner *cl = (UnitigCleaner*)arg;
  bool low_covg_snode = false, removable_tip = false;
  size_t i;

  CovgBuffer *cbuf = &cl->cbufs[threadid];
  fetch_coverages(nbuf, cbuf, cl->db_graph);

  // Covg is mean coverage of all kmers
  // size_t mean_covg, sum_covg = 0;
  // for(i = 0; i < cbuf->len; i++) sum_covg += cbuf->b[i];
  // mean_covg = sum_covg / cbuf->len;

  // Median coverage
  uint32_t median_covg = gca_median_uint32(cbuf->b, cbuf->len);

  low_covg_snode = (median_covg < cl->covg_threshold);

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
    // Keeping unitig
    for(i = 0; i < nbuf.len; i ++)
      (void)bitset_set_mt(cl->keep_flags, nbuf.b[i].key);

    // Update histograms
    update_kmer_covg_hist(cl->kmer_covgs_clean, cl->covg_arrsize,
                          cl->unitig_covg_clean, cl->covg_arrsize,
                          cl->len_hist_clean, cl->len_arrsize,
                          cbuf);
  }
}

/**
 * Remove unitigs with coverage < `covg_threshold` and tips shorter than
 * `min_keep_tip`.
 *
 * @param num_threads    Number of threads to use
 * @param covg_threshold Remove unitigs with mean covg < `covg_threshold`.
 *                       Ignored if 0.
 * @param min_keep_tip   Remove tips with length < `min_keep_tip`. Ignored if 0.
 * @param covgs_csv_path Path to write CSV of kmer coverage histogram
 * @param lens_csv_path  Path to write CSV of unitig length histogram
 *
 * `visited`, `keep` should each be at least db_graph.ht.capcity bits long
 *   and initialised to zero. On return,
 *   `visited` will be 1 at each original kmer index
 *   `keep` will be 1 at each retained kmer index
 **/
void clean_graph(size_t num_threads,
                 size_t covg_threshold, size_t min_keep_tip,
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

  if(covg_threshold > 0) {
    status("[cleaning] Removing unitigs with coverage < %zu...", covg_threshold);
    status("[cleaning]   Using kmer gamma method");
  }

  if(min_keep_tip > 0)
    status("[cleaning] Removing tips shorter than %zu...", min_keep_tip);

  status("[cleaning]   using %zu threads", num_threads);

  // Mark nodes to keep
  UnitigCleaner cl;
  unitig_cleaner_alloc(&cl, num_threads, covg_threshold,
                          min_keep_tip, keep, db_graph);
  supernodes_iterate(num_threads, visited, db_graph, unitig_mark, &cl);

  // Print numbers of kmers that are being removed

  char num_snodes_str[50], num_tips_str[50], num_tip_snodes_str[50];
  char num_snode_kmers_str[50], num_tip_kmers_str[50], num_tip_snode_kmers_str[50];
  ulong_to_str(cl.num_low_covg_snodes, num_snodes_str);
  ulong_to_str(cl.num_tips, num_tips_str);
  ulong_to_str(cl.num_tip_and_low_snodes, num_tip_snodes_str);
  ulong_to_str(cl.num_low_covg_snode_kmers, num_snode_kmers_str);
  ulong_to_str(cl.num_tip_kmers, num_tip_kmers_str);
  ulong_to_str(cl.num_tip_and_low_snode_kmers, num_tip_snode_kmers_str);

  status("[cleaning] Removing %s low coverage unitigs [%s kmer%s], "
         "%s unitig tips [%s kmer%s] "
         "and %s of both [%s kmer%s]",
         num_snodes_str,
         num_snode_kmers_str, util_plural_str(cl.num_low_covg_snode_kmers),
         num_tips_str,
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
    cleaning_write_covg_histogram(covgs_csv_path,
                                  cl.kmer_covgs_clean,
                                  cl.unitig_covg_clean,
                                  cl.covg_arrsize);
  }

  if(lens_csv_path != NULL) {
    cleaning_write_len_histogram(lens_csv_path,
                                 cl.len_hist_clean,
                                 cl.len_arrsize,
                                 db_graph->kmer_size);
  }

  unitig_cleaner_dealloc(&cl);
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
                                   const uint64_t *mean_covg_hist,
                                   size_t len)
{
  ctx_assert(len >= 2);
  ctx_assert(covg_hist[0] == 0);
  ctx_assert(mean_covg_hist[0] == 0);
  size_t i, end;

  FILE *fout = _open_histogram_file(path, "unitig coverage");
  if(fout == NULL) return;

  fprintf(fout, "Covg,NumKmers,NumUnitigs\n");
  for(end = len-1; end > 2 && covg_hist[end] == 0; end--) {}
  for(i = 1; i <= end; i++) {
    if(covg_hist[i] > 0)
      fprintf(fout, "%zu,%"PRIu64",%"PRIu64"\n", i, covg_hist[i], mean_covg_hist[i]);
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

  FILE *fout = _open_histogram_file(path, "unitig length");
  if(fout == NULL) return;

  fprintf(fout, "UnitigKmerLength,bp,Count\n");
  for(end = len-1; end > 1 && hist[end] == 0; end--) {}
  fprintf(fout, "1,%zu,%"PRIu64"\n", kmer_size, hist[1]);
  for(i = 2; i <= end; i++) {
    if(hist[i] > 0)
      fprintf(fout, "%zu,%zu,%"PRIu64"\n", i, kmer_size+i-1, hist[i]);
  }
  fclose(fout);
}
