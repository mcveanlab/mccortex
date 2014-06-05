#ifndef CLEAN_GRAPH_H_
#define CLEAN_GRAPH_H_

#include "db_graph.h"

// Get coverage threshold for removing supernodes
// If `min_keep_tip` is > 0, tips shorter than `min_keep_tip` are not used
// in measuring supernode coverage.
// `visited`, should each be at least db_graph.ht.capcity bits long
//   and initialised to zero
// `covgs_csv_path` and `lens_csv_path` are paths to files to write CSV
//   histogram of supernodes coverages and lengths BEFORE ANY CLEANING.
//   If NULL these are ignored.
Covg cleaning_get_threshold(size_t num_threads, double seq_depth,
                            const char *covgs_csv_path,
                            const char *lens_csv_path,
                            uint8_t *visited, dBGraph *db_graph);

// Remove low coverage supernodes and clip tips
// - Remove supernodes with coverage < `covg_threshold`
// - Remove tips shorter than `min_keep_tip`
// `visited`, `keep` should each be at least db_graph.ht.capcity bits long
//   and initialised to zero.
// `covgs_csv_path` and `lens_csv_path` are paths to files to write CSV
//   histogram of supernodes coverages and lengths AFTER CLEANING.
//   If NULL these are ignored.
void clean_graph(size_t num_threads, size_t covg_threshold, size_t min_keep_tip,
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
