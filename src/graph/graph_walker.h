#ifndef GRAPH_WALKER_H_
#define GRAPH_WALKER_H_

#include "db_graph.h"
#include "db_node.h"
#include "gpath_store.h"
#include "gpath_follow.h"

#include "madcrowlib/madcrow_list.h"

typedef struct
{
  uint8_t inleft:1, outleft:1, inright:1, outright:1;
  size_t num_nodes;
} GraphSection;

madcrow_list(gsec_list,GSecList,GraphSection);

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

// GraphWalker is not const because we update the path junction cache
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
void graph_walker_jump_along_snode(GraphWalker *wlk, hkey_t hkey,
                                   BinaryKmer bkmer, size_t num_nodes);

// return 1 on success, 0 otherwise
bool graph_walker_next(GraphWalker *wlk);
bool graph_walker_next_nodes(GraphWalker *wlk, size_t num_next,
                             const dBNode nodes[4], const Nucleotide bases[4]);

/**
 * Pick up counter paths for missing information check
 * @param prev_nodes nodes before wlk->node, oriented away from wlk->node
 *        i.e. the nodes you would reach if you walked from
 *        db_node_reverse(wlk->node)
*/
void graph_walker_add_counter_paths(GraphWalker *wlk,
                                    const dBNode prev_nodes[4],
                                    size_t num_prev);
// Prime for traversal
void graph_walker_prime(GraphWalker *wlk,
                        const dBNode *block, size_t n,
                        size_t max_context, bool forward,
                        size_t ctxcol, size_t ctpcol,
                        const dBGraph *db_graph);

bool graph_walker_agrees_contig(GraphWalker *wlk, const dBNode *block, size_t n,
                                bool forward);

#endif /* GRAPH_WALKER_H_ */
