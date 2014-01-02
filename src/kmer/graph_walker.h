#ifndef READ_PATH_H_
#define READ_PATH_H_

#include "path_store.h"
#include "db_graph.h"

typedef struct
{
  Nucleotide *bases;
  PathLen pos, len;
} FollowPath;


typedef struct
{
  // idx is -1 if failed, otherwise index of path
  int8_t idx;
  uint8_t status;
} GraphStep;

// GraphStep.status values:
#define GRPHWLK_FORWARD 0 /* Success: only one choice */
#define GRPHWLK_COLFWD 1 /* Success: only one choice in colour */
#define GRPHWLK_NOCOVG 2 /* Fail: no choices */
#define GRPHWLK_NOCOLCOVG 3 /* Fail: fork in pop but no choices in colour */
#define GRPHWLK_NOPATHS 4 /* Fail: fork in colour, no paths */
#define GRPHWLK_SPLIT_PATHS 5 /* Fail: fork in colour, paths split */
#define GRPHWLK_MISSING_PATHS 6 /* Fail: fork in colour, missing info */
#define GRPHWLK_USEPATH 7 /* Success: fork in colour, paths resolved */

// Was the last step resolving a split in this colour?
#define graphstep_is_fork(stp) ((stp).status > GRPHWLK_NOCOLCOVG)

typedef struct
{
  const dBGraph *const db_graph;
  const Colour ctxcol, ctpcol;

  // Current position
  hkey_t node;
  Orientation orient;
  BinaryKmer bkmer; // Oriented bkmer (i.e. not key)

  Nucleotide *data;
  FollowPath *allpaths;
  size_t max_path_len, max_num_paths;
  FollowPath **unused_paths, **curr_paths, **counter_paths;
  size_t num_unused, num_curr, num_new, num_counter;

  // Stats
  size_t fork_count;
  GraphStep last_step;
} GraphWalker;

// Need to pass number of colours in the graph
void graph_walker_alloc(GraphWalker *wlk);
void graph_walker_dealloc(GraphWalker *gw);

void graph_walker_print_state(const GraphWalker *wlk);

// Always call finish after calling init
void graph_walker_init(GraphWalker *wlk, const dBGraph *graph,
                       Colour ctxcol, Colour ctpcol,
                       hkey_t node, Orientation or);

void graph_walker_finish(GraphWalker *wlk);

// Hash a binary kmer + GraphWalker paths with offsets
uint32_t graph_walker_fasthash(const GraphWalker *wlk, const BinaryKmer bkmer);

// Returns index of choice or -1 along with status
GraphStep graph_walker_choose(const GraphWalker *wlk, size_t num_next,
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

void graph_walker_node_add_counter_paths(GraphWalker *wlk, Nucleotide prev_nuc);

#endif /* READ_PATH_H_ */
