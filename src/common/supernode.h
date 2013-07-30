#ifndef SUPERNODE_H_
#define SUPERNODE_H_

#include "graph_typedef.h"

typedef struct
{
  hkey_t start_node, end_node;
  Orientation start_orient, end_orient;
  uint32_t len;
} Supernode;

// Load a supernode from a given node
void supernode_load(hkey_t init_node, const dBGraph *db_graph, Supernode *supernode);

#endif
