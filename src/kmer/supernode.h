#ifndef SUPERNODE_H_
#define SUPERNODE_H_

#include "graph_typedef.h"

// Reverse order and orientations of nodes
void supernode_reverse(hkey_t *nlist, Orientation *olist, size_t len);

void supernode_normalise(hkey_t *nlist, Orientation *olist, size_t len);

// Extend a supernode, nlist[offset] and olist[offset] must already be set
// Walk along nodes starting from node/or, storing the supernode in nlist/olist
// Returns the number of nodes added, adds no more than `limit`
// return -1 if out of space and resize == false
int supernode_extend(const dBGraph *db_graph,
                     hkey_t **nlist, Orientation **olist,
                     size_t offset, boolean *cycle, size_t *arrlen,
                     boolean resize);

// Reallocates array if needs to resize
// returns length of supernode (always >=1)
// cycle is set to true if supernode is cycle
size_t supernode_find(dBGraph *db_graph, hkey_t node,
                      hkey_t **nlist, Orientation **olist,
                      boolean *cycle, size_t *arrlen);

void supernode_print(FILE *out, const dBGraph *db_graph,
                     hkey_t *nodes, Orientation *orients, size_t len);

void supernode_gzprint(gzFile out, const dBGraph *db_graph,
                       hkey_t *nodes, Orientation *orients, size_t len);

uint32_t supernode_read_starts(uint32_t *covgs, uint32_t len);

#endif
