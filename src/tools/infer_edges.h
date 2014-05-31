#ifndef INFER_EDGES_H_
#define INFER_EDGES_H_

// Return 1 if changed; 0 otherwise
bool infer_pop_edges(const BinaryKmer node_bkey, Edges *edges,
                     const Covg *covgs, const dBGraph *db_graph);

// Return 1 if changed; 0 otherwise
bool infer_all_edges(const BinaryKmer node_bkey, Edges *edges,
                     const Covg *covgs, const dBGraph *db_graph);

size_t infer_edges(size_t nthreads, bool add_all_edges, const dBGraph *db_graph);

#endif /* INFER_EDGES_H_ */
