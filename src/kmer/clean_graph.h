#ifndef CLEAN_GRAPH_H_
#define CLEAN_GRAPH_H_

#include "db_graph.h"

size_t cleaning_supernode_threshold(uint64_t *covgs, size_t len,
                                    double seq_depth, const dBGraph *db_graph);

// min_tip_len is the max tip length to KEEP
// visited should be zero'd before calling,
//   will have 1s for all remaining kmers on return
void cleaning_remove_tips(size_t min_tip_len, uint64_t *visited,
                          dBGraph *db_graph);

// visited should be zero'd before calling,
//   will have 1s for all remaining kmers on return
// Returns 0 if `covg_threshold` == 0 and computed threshold also <= 1
//  => no supernodes can be removed
Covg cleaning_remove_supernodes(bool clean, Covg covg_threshold,
                                double seq_depth, const char *dump_covgs,
                                uint64_t *visited, dBGraph *db_graph);

void cleaning_dump_covg_histogram(const char *path, uint64_t *hist, size_t len);

#endif /* CLEAN_GRAPH_H_ */
