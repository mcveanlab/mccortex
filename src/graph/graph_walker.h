#ifndef GRAPH_WALKER_H_
#define GRAPH_WALKER_H_

#include "db_graph.h"
#include "db_node.h"
#include "gpath_store.h"
#include "gpath_follow.h"

typedef struct
{
  uint8_t inleft:1, outleft:1, inright:1, outright:1;
  size_t num_nodes;
} GraphSection;

// Result from graph_walker_choose
typedef struct
{
  // idx is -1 if failed, otherwise index of node taken [0..3]
  int8_t idx;
  uint8_t status;
  bool node_has_col;
} GraphStep;

// GraphStep.status values:
#define GRPHWLK_FORWARD       0 /* Success: only one choice */
#define GRPHWLK_COLFWD        1 /* Success: only one choice in colour */
#define GRPHWLK_NOCOVG        2 /* Fail: no choices */
#define GRPHWLK_NOCOLCOVG     3 /* Fail: fork in pop but no choices in colour */
#define GRPHWLK_NOPATHS       4 /* Fail: fork in colour, no paths */
#define GRPHWLK_SPLIT_PATHS   5 /* Fail: fork in colour, paths split */
#define GRPHWLK_MISSING_PATHS 6 /* Fail: fork in colour, missing info */
#define GRPHWLK_USEPATH       7 /* Success: fork in colour, paths resolved */
#define GRPHWLK_NUM_STATES    8

#define GRPHWLK_FORWARD_STR       "GoForward"
#define GRPHWLK_COLFWD_STR        "GoColForward"
#define GRPHWLK_NOCOVG_STR        "FailNoCovg"
#define GRPHWLK_NOCOLCOVG_STR     "FailNoColCovg"
#define GRPHWLK_NOPATHS_STR       "FailNoPaths"
#define GRPHWLK_SPLIT_PATHS_STR   "FailSplitPaths"
#define GRPHWLK_MISSING_PATHS_STR "FailMissingPaths"
#define GRPHWLK_USEPATH_STR       "GoUsePath"

extern const char *graph_step_str[];

// Was the last step resolving a split in this colour?
#define graphstep_is_fork(stp) ((stp).status > GRPHWLK_NOCOLCOVG)

// Are we still walking?
#define grphwlk_status_is_good(stat) ((stat) <= GRPHWLK_COLFWD)

typedef struct
{
  const dBGraph *const db_graph;
  const GPathStore *const gpstore;
  const Colour ctxcol, ctpcol;

  // Current position
  dBNode node;
  BinaryKmer bkmer, bkey; // Oriented bkmer (i.e. not key) + hash bkey

  // Paths we are currently following
  GPathFollowBuffer paths, new_paths, cntr_paths;

  // Stats
  size_t fork_count;
  GraphStep last_step;
} GraphWalker;

// Get initial memory requirement
size_t graph_walker_est_mem();

// Need to pass number of colours in the graph
void graph_walker_alloc(GraphWalker *wlk);
void graph_walker_dealloc(GraphWalker *gw);

char* graph_walker_status2str(uint8_t status, char *str, size_t len);
void graph_walker_print_state_hist(const size_t arr[GRPHWLK_NUM_STATES]);
void graph_walker_print_state(const GraphWalker *wlk, FILE *fout);

// Always call finish after calling init
void graph_walker_init(GraphWalker *wlk, const dBGraph *graph,
                       Colour ctxcol, Colour ctpcol, dBNode node);

void graph_walker_finish(GraphWalker *wlk);

// Hash a binary kmer + GraphWalker paths with offsets
uint64_t graph_walker_hash64(GraphWalker *wlk);

// Returns index of choice or -1 along with status
GraphStep graph_walker_choose(GraphWalker *wlk, size_t num_next,
                              const dBNode next_nodes[4],
                              const Nucleotide next_bases[4]);

// Move to the next node
// If fork is true, node is the result of taking a fork -> slim down paths
void graph_walker_force(GraphWalker *wlk, hkey_t hkey, Nucleotide base,
                        bool fork);

// Jump to a new node within the current sample supernode
// (can actually be any node up until the end of the current supernode)
void graph_walker_jump_along_snode(GraphWalker *wlk, hkey_t hkey, BinaryKmer bkmer);

// return 1 on success, 0 otherwise
bool graph_walker_next(GraphWalker *wlk);
bool graph_walker_next_nodes(GraphWalker *wlk, size_t num_next,
                          const dBNode nodes[4], const Nucleotide bases[4]);

void graph_walker_add_counter_paths(GraphWalker *wlk,
                                    hkey_t prev_nodes[4],
                                    Orientation prev_orients[4],
                                    size_t num_prev);


// Fast traversal of a list of nodes using the supplied GraphWalker
// Only visits nodes deemed informative + last node
// Must have previously initialised or walked to the prior node,
// using: graph_walker_init, graph_walker_force, graph_walker_jump_along_snode,
// graph_traverse or graph_walker_next_nodes
// i.e. wlk->node is a node adjacent to arr[0]
void graph_walker_fast_traverse(GraphWalker *wlk, const dBNode *nodes, size_t n,
                                bool forward);

// Force traversal of every node
void graph_walker_slow_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                bool forward);

// Prime for traversal
void graph_walker_prime(GraphWalker *wlk,
                        const dBNode *block, size_t n,
                        size_t max_context, bool forward,
                        size_t ctxcol, size_t ctpcol,
                        const dBGraph *db_graph);

bool graph_walker_agrees_contig(GraphWalker *wlk, const dBNode *block, size_t n,
                                bool forward);

// How many junctions are left to be traversed in our longest remaining path
size_t graph_walker_get_max_path_junctions(const GraphWalker *wlk);

#endif /* GRAPH_WALKER_H_ */
