#ifndef CLEAN_GRAPH_H_
#define CLEAN_GRAPH_H_

#include "db_graph.h"

/**
 * Pick a cleaning threshold from kmer coverage histogram. Assumes low coverage
 * kmers are all due to error, to which it fits a gamma distribution. Then
 * chooses a cleaning threshold such than FDR (uncleaned kmers) occur at a rate
 * of < the FDR paramater.
 *
 * Translated from Gil McVean's proposed method in R code
 *
 * @param kmer_covg Histogram of kmer counts at coverages 1,2,.. arrlen-1
 * @param arrlen    Length of array kmer_covg
 * @param fdr_limit False discovery rate for a single kmer coverage
 *                  (1/1000 i.e. 0.001 is reasonable)
 * @param alpha_est_ptr If not NULL, used to return estimate for alpha
 * @param beta_est_ptr  If not NULL, used to return estimate for beta
 * @return -1 if no cut-off satisfies FDR, otherwise returs coverage cutoff
 */
int cleaning_pick_kmer_threshold(const uint64_t *kmer_covg, size_t arrlen,
                                 double fdr_limit,
                                 double *alpha_est_ptr, double *beta_est_ptr);

/**
 * Calculate cleaning threshold for supernodes from a given distribution
 * of supernode coverages
 * @param covgs histogram of supernode coverages
 */
size_t cleaning_pick_supernode_threshold(const uint64_t *covgs, size_t len,
                                         double seq_depth,
                                         const dBGraph *db_graph);

/**
 * Get coverage threshold for removing supernodes
 *
 * @param visited should be at least db_graph.ht.capcity bits long and initialised
 *                to zero. On return, it will be 1 at each original kmer index
 * @param covgs_csv_path
 * @param lens_csv_path  paths to files to write CSV histogram of supernodes
                         coverages and lengths BEFORE ANY CLEANING.
 *                       If NULL these are ignored.
 * @return threshold to clean or -1 on error
 */
int cleaning_get_threshold(size_t num_threads,
                           const char *covgs_csv_path,
                           const char *lens_csv_path,
                           uint8_t *visited,
                           const dBGraph *db_graph);

/**
 * Remove low coverage supernodes and clip tips
 * - Remove supernodes with coverage < `covg_threshold`
 * - Remove tips shorter than `min_keep_tip`
 * `visited`, `keep` should each be at least db_graph.ht.capcity bits long
 *   and initialised to zero.
 * `covgs_csv_path` and `lens_csv_path` are paths to files to write CSV
 *   histogram of supernodes coverages and lengths AFTER CLEANING.
 *   If NULL these are ignored.
 */
void clean_graph(size_t num_threads,
                 size_t covg_threshold, size_t min_keep_tip,
                 const char *covgs_csv_path, const char *lens_csv_path,
                 uint8_t *visited, uint8_t *keep, dBGraph *db_graph);

void cleaning_write_covg_histogram(const char *path,
                                   const uint64_t *covg_hist,
                                   const uint64_t *kmer_hist,
                                   size_t len);

void cleaning_write_len_histogram(const char *path,
                                  const uint64_t *hist, size_t len,
                                  size_t kmer_size);

#endif /* CLEAN_GRAPH_H_ */
