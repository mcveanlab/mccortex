#ifndef READ_PATH_H_
#define READ_PATH_H_

#include "binary_paths.h"
#include "db_graph.h"

typedef struct
{
  Nucleotide *bases;
  PathLen pos, len;
} FollowPath;

typedef struct
{
  const dBGraph *const db_graph;
  const Colour colour;

  // Current position
  hkey_t node;
  Orientation orient;
  BinaryKmer bkmer; // Oriented bkmer (i.e. not key)

  Nucleotide *data;
  FollowPath *allpaths;
  size_t max_path_len, max_num_paths;
  FollowPath **unused_paths, **curr_paths, **counter_paths;
  size_t num_unused, num_curr, num_new, num_counter;
} GraphWalker;

// Need to pass number of colours in the graph
void graph_walker_alloc(GraphWalker *wlk);
void graph_walker_dealloc(GraphWalker *gw);

// Always call finish after calling init
void graph_walker_init(GraphWalker *wlk, const dBGraph *graph, Colour colour,
                       hkey_t node, Orientation or);

// #define MAX_WALK_BACK_NODES 100

// // context is now many nodes to go back (up to MAX_WALK_BACK_NODES)
// // Remember to call finish when done with wlk
// void graph_walker_init_context(GraphWalker *wlk, const dBGraph *db_graph,
//                                uint64_t *visited, Colour colour,
//                                hkey_t node, Orientation orient);

void graph_walker_finish(GraphWalker *wlk);

// Returns index of choice or -1
int graph_walker_choose(const GraphWalker *wlk, size_t num_next,
                        const hkey_t next_nodes[4],
                        const Nucleotide next_bases[4]);

// Move to the next node
// If fork is true, node is the result of taking a fork -> slim down paths
void graph_traverse_force(GraphWalker *wlk, hkey_t node, Nucleotide base,
                          boolean fork);

// Jump to a new node (any node up until the end of the current supernode or the
// first node of the next supernode)
void graph_traverse_force_jump(GraphWalker *wlk, hkey_t node, BinaryKmer bkmer,
                               boolean fork);

// return 1 on success, 0 otherwise
boolean graph_traverse(GraphWalker *wlk);
boolean graph_traverse_nodes(GraphWalker *wlk, size_t num_next,
                             const hkey_t nodes[4], const Nucleotide bases[4]);

void graph_walker_add_counter_paths(GraphWalker *wlk,
                                    hkey_t prev_nodes[4],
                                    Orientation prev_orients[4],
                                    size_t num_prev);

void graph_walker_node_add_counter_paths(GraphWalker *wlk,
                                         hkey_t node, Orientation orient,
                                         Nucleotide prev_nuc);

#endif /* READ_PATH_H_ */
