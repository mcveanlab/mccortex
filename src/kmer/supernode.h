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

// Reallocates array if needs to resize
// returns length of supernode (always >=1)
void supernode_find(hkey_t node, dBNodeBuffer *nbuf, const dBGraph *db_graph);

// Count number of read starts using coverage data
uint32_t supernode_read_starts(const uint32_t *covgs, uint32_t len);

void supernode_write_len_distrib(FILE *fout, const char *path, size_t histlen,
                                 uint64_t *visited, const dBGraph *db_graph);

#endif
