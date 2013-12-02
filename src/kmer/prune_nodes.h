#ifndef PRUNE_NODES_H_
#define PRUNE_NODES_H_

//
// Pruning nodes from the graph
//
#include "graph_typedef.h"
#include "db_node.h"

void prune_node(dBGraph *db_graph, hkey_t node);

// Supernode pruning used by ctx_clean
void prune_supernode(dBGraph *db_graph, dBNode *nodes, size_t len);

// Used by ctx_subgraph
// flags is a bit array, one bit per kmer
void prune_nodes_lacking_flag(dBGraph *graph, uint64_t *flags);

// Currently unused
// remove nodes if not in any colour
// i.e. db_node_has_col(graph,node,colour) == false for all colours
void prune_uncoloured_nodes(dBGraph *db_graph);

#endif /* PRUNE_NODES_H_ */
