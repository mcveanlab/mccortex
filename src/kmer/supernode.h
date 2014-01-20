#ifndef SUPERNODE_H_
#define SUPERNODE_H_

#include "cortex_types.h"
#include "db_node.h"

// Reverse order and orientations of nodes
void supernode_reverse(dBNode *nlist, size_t len);

void supernode_normalise(dBNode *nlist, size_t len);

// Extend a supernode, nlist[offset] and olist[offset] must already be set
// Walk along nodes starting from node/or, storing the supernode in nlist/olist
// Returns the number of nodes added, adds no more than `limit`
// return -1 if out of space and resize == false
int supernode_extend(dBNode **nlist, size_t offset, size_t *arrlen,
                     boolean resize, const dBGraph *db_graph);

// Reallocates array if needs to resize
// returns length of supernode (always >=1)
size_t supernode_find(hkey_t node, dBNode **nlist, size_t *arrlen,
                      const dBGraph *db_graph);

uint32_t supernode_read_starts(uint32_t *covgs, uint32_t len);

#endif
