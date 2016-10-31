#ifndef INFER_EDGES_H_
#define INFER_EDGES_H_

#include "cortex_types.h"
#include "db_graph.h"

// `pop_edges` if true, only add edges that are in at least one other colour
//  -> If two kmers are in a sample and the population has an edges between
//     them, add edge to sample.
// Return 1 if changed; 0 otherwise
bool infer_kmer_edges(const BinaryKmer node_bkey, bool pop_edges,
                      Edges *edges, const Covg *covgs,
                      const dBGraph *db_graph);

size_t infer_edges(size_t nthreads, bool add_all_edges, const dBGraph *db_graph);

#endif /* INFER_EDGES_H_ */
