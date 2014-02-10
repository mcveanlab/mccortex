#ifndef GRAPH_WALKER_H_
#define GRAPH_WALKER_H_

#include "db_graph.h"
#include "db_node.h"
#include "path_store.h"

// This struct is packed so we can hash it quickly
struct FollowPathStruct
{
  const uint8_t *seq;
  PathLen pos, len;
  // A small buffer of upcoming 24 bases
  PathLen first_cached; // first base in buffer (multiple of 4: 0,4,8,...)
  uint8_t cache[6]; // first..first+24-1 (24 bases)
} __attribute__((packed));

typedef struct FollowPathStruct FollowPath;

FollowPath follow_path_create(const uint8_t *seq, PathLen plen);

#include "objbuf_macro.h"
create_objbuf(path_buf,PathBuffer,FollowPath)

// Result from graph_walker_choose
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
  const PathStore *const pstore;
  const Colour ctxcol, ctpcol;

  // Current position
  // hkey_t node;
  // Orientation orient;
  dBNode node;
  BinaryKmer bkmer, bkey; // Oriented bkmer (i.e. not key) + hash bkey

  // Paths we are currently following
  PathBuffer paths, new_paths, cntr_paths;

  // Stats
  size_t fork_count;
  GraphStep last_step;
} GraphWalker;

// Get initial memory requirement
size_t graph_walker_est_mem();

// Need to pass number of colours in the graph
void graph_walker_alloc(GraphWalker *wlk);
void graph_walker_dealloc(GraphWalker *gw);

void graph_walker_print_state(const GraphWalker *wlk);

// Always call finish after calling init
void graph_walker_init(GraphWalker *wlk, const dBGraph *graph,
                       Colour ctxcol, Colour ctpcol, dBNode node);

void graph_walker_finish(GraphWalker *wlk);

// Hash a binary kmer + GraphWalker paths with offsets
// uint32_t graph_walker_hash(const GraphWalker *wlk);
uint64_t graph_walker_hash64(const GraphWalker *wlk);

// Returns index of choice or -1 along with status
GraphStep graph_walker_choose(const GraphWalker *wlk, size_t num_next,
                              const dBNode next_nodes[4],
                              const Nucleotide next_bases[4]);

// Move to the next node
// If fork is true, node is the result of taking a fork -> slim down paths
void graph_traverse_force(GraphWalker *wlk, hkey_t hkey, Nucleotide base,
                          boolean fork);

// Jump to a new node within the current sample supernode
// (can actually be any node up until the end of the current supernode)
void graph_walker_jump_snode_end(GraphWalker *wlk, hkey_t hkey, BinaryKmer bkmer);

// return 1 on success, 0 otherwise
boolean graph_traverse(GraphWalker *wlk);
boolean graph_traverse_nodes(GraphWalker *wlk, size_t num_next,
                             const dBNode nodes[4], const Nucleotide bases[4]);

void graph_walker_add_counter_paths(GraphWalker *wlk,
                                    hkey_t prev_nodes[4],
                                    Orientation prev_orients[4],
                                    size_t num_prev);


// Fast traversal of a list of nodes using the supplied GraphWalker
// Only visits nodes deemed informative + last node
// Must have previously initialised or walked to the prior node,
// using: graph_walker_init, graph_traverse_force, graph_walker_jump_snode_end,
// graph_traverse or graph_traverse_nodes
// i.e. wlk->node is a node adjacent to arr[0]
void graph_walker_fast_traverse(GraphWalker *wlk, const dBNode *nodes, size_t n,
                                boolean forward);

// Force traversal of every node
void graph_walker_slow_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                boolean forward);

// Prime for traversal
void graph_walker_prime(GraphWalker *wlk,
                        const dBNode *block, size_t n,
                        size_t max_context, boolean forward,
                        size_t ctxcol, size_t ctpcol,
                        const dBGraph *db_graph);

#endif /* GRAPH_WALKER_H_ */
