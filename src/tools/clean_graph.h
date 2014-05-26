#ifndef CLEAN_GRAPH_H_
#define CLEAN_GRAPH_H_

#include "db_graph.h"

// Get coverage threshold for removing supernodes
// If `min_keep_tip` is > 0, tips shorter than `min_keep_tip` are not used
// in measuring supernode coverage.
// `visited`, should each be at least db_graph.ht.capcity bits long
//   and initialised to zero
Covg cleaning_get_threshold(size_t num_threads, size_t min_keep_tip,
                            double seq_depth, const char *dump_covgs,
                            uint8_t *visited, dBGraph *db_graph);

// Remove low coverage supernodes and clip tips
// - Remove supernodes with coverage < `covg_threshold`
// - Remove tips shorter than `min_keep_tip`
// `visited`, `keep` should each be at least db_graph.ht.capcity bits long
//   and initialised to zero.
void clean_graph(size_t num_threads, size_t covg_threshold, size_t min_keep_tip,
                 uint8_t *visited, uint8_t *keep, dBGraph *db_graph);

void cleaning_dump_covg_histogram(const char *path, uint64_t *hist, size_t len);

#endif /* CLEAN_GRAPH_H_ */
