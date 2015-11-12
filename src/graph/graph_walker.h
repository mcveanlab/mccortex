#ifndef GRAPH_WALKER_H_
#define GRAPH_WALKER_H_

#include "db_graph.h"
#include "db_node.h"
#include "gpath_store.h"
#include "gpath_follow.h"
#include "graph_step.h"

#include "madcrowlib/madcrow_list.h"

typedef struct
{
  uint32_t num_nodes:30, in_fork:1, out_fork:1;
} GraphSegment;

madcrow_list(gseg_list,GSegList,GraphSegment);

typedef struct
{
  const dBGraph *db_graph;
  const GPathStore *gpstore;
  Colour ctxcol, ctpcol;
  bool missing_path_check; // if true do missing path check
  size_t *used_paths; // bit array for used paths; cast volatile when used
  // 1 bit for each path loaded

  // Current position
  dBNode node;
  BinaryKmer bkey;

  // Paths we are currently following
  GPathFollowBuffer paths, cntr_paths;
  GSegList gsegs;

  // Statistics
  size_t fork_count; // how many forks we have traversed
  GraphStep last_step;
} GraphWalker;

void graph_walker_print_state(const GraphWalker *wlk, FILE *fout);

// Get initial memory requirement
size_t graph_walker_est_mem();

// Need to pass number of colours in the graph
void graph_walker_alloc(GraphWalker *wlk, const dBGraph *graph);
void graph_walker_dealloc(GraphWalker *gw);

void graph_walker_setup(GraphWalker *wlk, bool missing_path_check,
                        Colour ctxcol, Colour ctpcol,
                        const dBGraph *graph);

// Always call finish after calling start
void graph_walker_start(GraphWalker *wlk, dBNode node);
void graph_walker_finish(GraphWalker *wlk);

// Hash a binary kmer + GraphWalker paths with offsets
uint64_t graph_walker_hash64(GraphWalker *wlk);

/**
 * Make a choice at a junction
 * @return index of choice or -1
 */
GraphStep graph_walker_choose(const GraphWalker *wlk, size_t num_next,
                              const dBNode next_nodes[4],
                              const Nucleotide next_bases[4]);

/**
 * Move to the next node
 * @param is_fork If true, node is the result of taking a fork (updates paths)
 */
void graph_walker_force(GraphWalker *wlk, dBNode node, bool is_fork);

/**
 * Jump to a new node within the current sample unitig
 * (can actually be any node up until the end of the current unitig)
 * @param num_nodes is number of nodes we have moved forward
 */
void graph_walker_jump_along_snode(GraphWalker *wlk, dBNode node, size_t num_nodes);

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

/**
 * Traversal of every node in a list of nodes using the supplied GraphWalker
 * Visits each node specifed
 */
void graph_walker_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                           bool forward);

// Prime for traversal
void graph_walker_prime(GraphWalker *wlk,
                        const dBNode *block, size_t n,
                        size_t max_context, bool forward);

/**
 * Check the graph walker doesn't veer away from the given contig
 * At each node:
 *  a. If we can't progress -> success
 *  b. If we can and it doesn't match what we expected -> disagrees
 * @param forward Traverse contig forward, otherwise reverse complement nodes
 *                 and work backwards
 */
bool graph_walker_agrees_contig(GraphWalker *wlk, const dBNode *block, size_t n,
                                bool forward);

#endif /* GRAPH_WALKER_H_ */
