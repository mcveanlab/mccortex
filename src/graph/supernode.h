#ifndef SUPERNODE_H_
#define SUPERNODE_H_

#include "cortex_types.h"
#include "db_graph.h"
#include "db_node.h"

// Reverse order and orientations of nodes
void db_nodes_reverse_complement(dBNode *nlist, size_t len);

// Orient supernode
// Once oriented, supernode has lowest poosible kmerkey at the beginning,
// oriented FORWARDs if possible
void supernode_normalise(dBNode *nlist, size_t len, const dBGraph *db_graph);

// Extend a supernode, nlist[offset] and olist[offset] must already be set
// Walk along nodes starting from node/or, storing the supernode in nlist/olist
// Returns the number of nodes added, adds no more than `limit`
// return false if out of space and limit > 0
bool supernode_extend(dBNodeBuffer *nbuf, size_t limit,
                      const dBGraph *db_graph);

// Fills with supernode that contains hkey
// Does not reset nbuf
void supernode_find(hkey_t node, dBNodeBuffer *nbuf, const dBGraph *db_graph);

// Count number of read starts using coverage data
size_t supernode_read_starts(const Covg *covgs, size_t len);
size_t supernode_covg_mean(const Covg *covgs, size_t len);

void supernodes_iterate(size_t nthreads, uint8_t *visited,
                        const dBGraph *db_graph,
                        void (*func)(dBNodeBuffer _nbuf,
                                     size_t threadid,
                                     void *_arg),
                        void *arg);

#endif /* SUPERNODE_H_ */
