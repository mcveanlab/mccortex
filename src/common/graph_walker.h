#ifndef READ_PATH_H_
#define READ_PATH_H_

#include "binary_paths.h"
#include "db_graph.h"

typedef struct
{
  dBGraph *const db_graph;
  const Colour colour; // This is not currently used

  // Current position
  hkey_t node;
  BinaryKmer bkmer;
  Orientation orient;

  // Current paths
  path_t **curr_paths, *paths_data;
  size_t num_paths, paths_capacity;

  uint64_t *prev_paths;
  size_t num_pp, pp_cap;
} GraphWalker;

void graph_walker_alloc(GraphWalker *wlk);
void graph_walker_dealloc(GraphWalker *gw);

void graph_walker_init(GraphWalker *wlk, dBGraph *graph, Colour colour,
                       hkey_t node, Orientation or);

void graph_walker_finish(GraphWalker *wlk);

// DEV: write this function
// void graph_walker_init_context(GraphWalker *wlk, dBGraph *graph, int colour,
//                                Element **els, Orientation *ors, int len);

// Index of choice or -1
int graph_walker_choose(GraphWalker *wlk, size_t num_next, Nucleotide bases[4]);

// If fork is true, node is the result of taking a fork -> slim down paths
void graph_traverse_force(GraphWalker *wlk, hkey_t node, BinaryKmer bkmer,
                          Orientation orient, boolean fork);

// return 1 on success, 0 otherwise
boolean graph_traverse(GraphWalker *wlk);
boolean graph_traverse_nodes(GraphWalker *wlk, size_t num_next,
                             hkey_t nodes[4], BinaryKmer bkmers[4],
                             Orientation orients[4]);

#endif /* READ_PATH_H_ */
