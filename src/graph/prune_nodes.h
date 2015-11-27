#ifndef PRUNE_NODES_H_
#define PRUNE_NODES_H_

//
// Pruning nodes from the graph
//
#include "cortex_types.h"
#include "db_node.h"

// Remove a node from the graph, do not edit any edges / adjacent nodes
// Threadsafe
void prune_node_without_edges_mt(dBGraph *db_graph, hkey_t hkey);

void prune_node(dBGraph *db_graph, hkey_t node);

// Unitig pruning used by ctx_clean
void prune_unitig(dBNode *nodes, size_t len, dBGraph *db_graph);

// Used by ctx_subgraph.c, clean_graph.c
// flags is a bit array, one bit per kmer
void prune_nodes_lacking_flag(size_t num_threads, const uint8_t *flags,
                              dBGraph *db_graph);

// Currently unused
// remove nodes if not in any colour
// i.e. db_node_has_col(graph,node,colour) == false for all colours
void prune_uncoloured_nodes(dBGraph *db_graph);

#endif /* PRUNE_NODES_H_ */
