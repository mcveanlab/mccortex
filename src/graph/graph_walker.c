#include "global.h"
#include "util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"
#include "gpath.h"
#include "binary_seq.h"

// hash functions
#include "misc/jenkins.h"
#include "misc/twang.h"
#include "misc/city.h"

#ifdef CTXVERBOSE
#define DEBUG_WALKER 1
#endif

// How many junctions are left to be traversed in our longest remaining path
static size_t graph_walker_get_max_path_junctions(const GraphWalker *wlk)
{
  size_t i, rem, max = 0;
  for(i = 0; i < wlk->paths.len; i++) {
    rem = (size_t)(wlk->paths.b[i].len - wlk->paths.b[i].pos);
    max = MAX2(rem, max);
  }
  return max;
}

static void print_path_list(const GPathFollowBuffer *pbuf, FILE *fout)
{
  size_t i, j;
  GPathFollow *path;

  for(i = 0; i < pbuf->len; i++) {
    path = &pbuf->b[i];
    fprintf(fout, "   %p ", path->gpath->seq);
    for(j = 0; j < path->len; j++)
      fputc(dna_nuc_to_char(gpath_follow_get_base(path, j)), fout);
    fprintf(fout, " [%zu/%zu] age: %zu %c\n", (size_t)path->pos,
            (size_t)path->len, (size_t)path->age,
            gpath_follow_get_base(path, path->pos));
  }
}

void graph_walker_print_state(const GraphWalker *wlk, FILE *fout)
{
  char bkeystr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, bkeystr);
  fprintf(fout, " GWState: %s:%i ctx: %zu ctp: %zu\n",
         bkeystr, wlk->node.orient, wlk->ctxcol, wlk->ctpcol);
  fprintf(fout, "  num_curr: %zu\n", wlk->paths.len);
  print_path_list(&wlk->paths, fout);
  fprintf(fout, "  num_counter: %zu\n", wlk->cntr_paths.len);
  print_path_list(&wlk->cntr_paths, fout);
  fprintf(fout, "--\n");
}

size_t graph_walker_est_mem()
{
  return sizeof(GPathFollow)*1024;
}

// Allocate memory, default to colour 0 and using missing info check
void graph_walker_alloc(GraphWalker *wlk, const dBGraph *graph)
{
  memset(wlk, 0, sizeof(GraphWalker));
  gpath_follow_buf_alloc(&wlk->paths, 256);
  gpath_follow_buf_alloc(&wlk->cntr_paths, 512);
  gseg_list_alloc(&wlk->gsegs, 128);
  graph_walker_setup(wlk, true, 0, 0, graph);
}

// Free memory
void graph_walker_dealloc(GraphWalker *wlk)
{
  gpath_follow_buf_dealloc(&wlk->paths);
  gpath_follow_buf_dealloc(&wlk->cntr_paths);
  gseg_list_dealloc(&wlk->gsegs);
  memset(wlk, 0, sizeof(GraphWalker));
}

void graph_walker_setup(GraphWalker *wlk, bool missing_path_check,
                        Colour ctxcol, Colour ctpcol,
                        const dBGraph *graph)
{
  // Check that the graph is loaded properly (all edges merged into one colour)
  ctx_assert(graph->num_edge_cols == 1);
  ctx_assert(graph->num_of_cols == 1 || graph->node_in_cols != NULL);

  wlk->db_graph = graph;
  wlk->gpstore = &graph->gpstore;
  wlk->ctxcol = ctxcol;
  wlk->ctpcol = ctpcol;
  wlk->missing_path_check = missing_path_check;
}

static inline void _gw_gseg_init(GraphWalker *wlk)
{
  ctx_assert(gseg_list_len(&wlk->gsegs) == 0);
  GraphSegment gseg = {.in_fork = false, .out_fork = false, .num_nodes = 1};
  gseg_list_unshift(&wlk->gsegs, &gseg, 1);
}

static inline void _gw_gseg_update(GraphWalker *wlk,
                                   bool fw_fork, bool rv_fork,
                                   size_t num_nodes)
{
  // First GraphSection is the one the most recent one we saw
  size_t i, len = gseg_list_len(&wlk->gsegs);
  GraphSegment *first = gseg_list_getptr(&wlk->gsegs, 0);
  ctx_assert(len > 0);

  // Previous section may have ended in a fork - update it
  first->out_fork |= fw_fork;

  if(fw_fork || rv_fork)
  {
    // Start a new section
    ctx_assert(num_nodes == 1);
    GraphSegment gseg = {.in_fork = rv_fork, .out_fork = 0, .num_nodes = 0};
    gseg_list_unshift(&wlk->gsegs, &gseg, 1);

    len = gseg_list_len(&wlk->gsegs);
    first = gseg_list_getptr(&wlk->gsegs, 0);

    // Update all path ages: age is the graph section index
    for(i = 0; i < wlk->paths.len; i++)
      wlk->paths.b[i].age++;

    for(i = 0; i < wlk->cntr_paths.len; i++)
      wlk->cntr_paths.b[i].age++;

    // Drop graph sections without any paths (apart from most recent one)
    size_t max_segs = 1;

    if(wlk->paths.len)
      max_segs = MAX2(max_segs, wlk->paths.b[0].age + 1);
    if(wlk->cntr_paths.len)
      max_segs = MAX2(max_segs, wlk->cntr_paths.b[0].age + 1);

    ctx_assert2(max_segs <= len, "%zu > %zu", max_segs, len);
    gseg_list_pop(&wlk->gsegs, NULL, len - max_segs);
  }

  // printf("jump %zu\n", num_nodes);
  first->num_nodes += num_nodes;
}

// Returns number of paths picked up
// next_nuc only used if counter == true and node has out-degree > 1
static inline size_t pickup_paths(GraphWalker *wlk, dBNode node,
                                  bool counter, Nucleotide next_nuc)
{
  const dBGraph *db_graph = wlk->db_graph;
  const GPathStore *gpstore = wlk->gpstore;
  GPathFollowBuffer *pbuf = counter ? &wlk->cntr_paths : &wlk->paths;
  const size_t ncols = gpstore->gpset.ncols, num_paths = pbuf->len;

  // Picking up paths is turned off
  if(!gpath_store_use_traverse(gpstore)) return 0;
  if(!db_node_in_col(db_graph, wlk->node.key, wlk->ctxcol)) return 0;

  // DEBUG
  // char kstr[MAX_KMER_SIZE+3]; // <kmer>:<orient>
  // db_node_to_str(wlk->db_graph, node, kstr);
  // printf("pickup %s paths from: %s\n", counter ? "cntr" : "curr", kstr);
  // END DEBUG

  // cntr_filter_nuc0 is needed for picking up counter paths with outdegree > 1
  bool cntr_filter_nuc0
    = (counter && db_node_outdegree_in_col(node, wlk->ctxcol, db_graph) > 1);

  GPath *gpath = gpath_store_fetch_traverse(gpstore, node.key);

  for(; gpath != NULL; gpath = gpath->next)
  {
    if(node.orient == gpath->orient && gpath_has_colour(gpath, ncols, wlk->ctpcol))
    {
      GPathFollow fpath = gpath_follow_create(gpath);

      if(!cntr_filter_nuc0) gpath_follow_buf_add(pbuf, fpath);
      else if(gpath_follow_get_base(&fpath, 0) == next_nuc) {
        // Loading a counter path at a fork
        fpath.pos++; // already took a base
        // check there are still junctions to take
        if(fpath.pos < fpath.len)
          gpath_follow_buf_add(pbuf, fpath);
      }
    }
  }

  #ifdef DEBUG_WALKER
    char bkey_str[MAX_KMER_SIZE+1], node_str[MAX_KMER_SIZE+1];
    BinaryKmer node_bkey = db_node_get_bkey(wlk->db_graph, node.key);
    binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, bkey_str);
    binary_kmer_to_str(node_bkey, wlk->db_graph->kmer_size, node_str);
    status("  pickup_paths(): %s:%i node:%s:%i picked up %zu %s paths cntr_filter_nuc0:%i",
           bkey_str, wlk->node.orient, node_str, node.orient,
           pbuf->len - num_paths, counter ? "counter" : "forward",
           cntr_filter_nuc0);
  #endif

  // DEBUG
  // if(pbuf->len > num_paths) {
  //   fprintf(stderr, "  Picked up %zu paths; node: %zu %s\n",
  //           pbuf->len - num_paths, (size_t)node.key, counter ? " [cntr]" : "");
  // }

  return pbuf->len - num_paths;
}

/**
 * Pick up counter paths for missing information check
 * @param prev_nodes nodes before wlk->node, oriented towards wlk->node
*/
void graph_walker_add_counter_paths(GraphWalker *wlk,
                                    const dBNode prev_nodes[4],
                                    size_t num_prev)
{
  if(!wlk->missing_path_check) return;

  size_t i;
  Nucleotide next_base = bkmer_get_last_nuc(wlk->bkey, wlk->node.orient,
                                            wlk->db_graph->kmer_size);

  // Reverse orientation, pick up paths
  for(i = 0; i < num_prev; i++)
    pickup_paths(wlk, prev_nodes[i], true, next_base);
}


void graph_walker_start(GraphWalker *wlk, dBNode node)
{
  ctx_assert(wlk->paths.len == 0);
  ctx_assert(wlk->cntr_paths.len == 0);
  ctx_assert(gseg_list_len(&wlk->gsegs) == 0);

  wlk->node = node;

  // stats
  wlk->fork_count = 0;
  memset(&wlk->last_step, 0, sizeof(wlk->last_step));
  wlk->last_step.idx = -1;

  // Get binary kmer
  wlk->bkey = db_node_get_bkey(wlk->db_graph, node.key);

  #ifdef DEBUG_WALKER
    char kmer_str[MAX_KMER_SIZE+1];
    binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, kmer_str);
    status("  graph_walker_start(): %s:%i", kmer_str, wlk->node.orient);
  #endif

  _gw_gseg_init(wlk);

  // Pick up new paths
  pickup_paths(wlk, wlk->node, false, 0);

  ctx_assert(gseg_list_len(&wlk->gsegs) == 1);
  // printf("zero: %u\n", gseg_list_getptr(&wlk->gsegs, 0)->num_nodes);

  // Don't pick up counter paths on init() - there is no point
  // since there is no scenario where they'd help if picked up here
  // Edges edges = db_node_get_edges(db_graph, wlk->node.key, 0);
  // if(edges_get_indegree(edges, wlk->node.orient) > 1)
  //   _graph_walker_pickup_counter_paths(wlk, lost_nuc);
}

void graph_walker_finish(GraphWalker *wlk)
{
  #ifdef DEBUG_WALKER
    char kmer_str[MAX_KMER_SIZE+1];
    binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, kmer_str);
    status("  graph_walker_finish(): %s:%i", kmer_str, wlk->node.orient);
  #endif

  gpath_follow_buf_reset(&wlk->paths);
  gpath_follow_buf_reset(&wlk->cntr_paths);
  gseg_list_reset(&wlk->gsegs);
}

//
// Hash function
//

uint64_t graph_walker_hash64(GraphWalker *wlk)
{
  size_t path_bytes, cntr_path_bytes;
  path_bytes = wlk->paths.len*sizeof(GPathFollow);
  cntr_path_bytes = wlk->cntr_paths.len*sizeof(GPathFollow);

  uint64_t hash = ((uint64_t)wlk->node.key<<1) | wlk->node.orient;
  hash = CityHash64WithSeeds((const char*)wlk->paths.b, path_bytes,
                             hash, wlk->paths.len);
  hash = CityHash64WithSeeds((const char*)wlk->cntr_paths.b, cntr_path_bytes,
                             hash, wlk->cntr_paths.len);

  return hash;
}

//
// Junction decision
//

static inline void update_path_forks(const GPathFollowBuffer *pbuf, bool taken[4])
{
  size_t i;
  GPathFollow *path;
  for(i = 0; i < pbuf->len; i++) {
    path = &pbuf->b[i];
    taken[gpath_follow_get_base(path, path->pos)] = true;
  }
}

static inline void _corrupt_paths(const GraphWalker *wlk, size_t num_next,
                                  const dBNode nodes[4],
                                  const Nucleotide bases[4])
__attribute__((noreturn));

// If we call this funcitons, something has gone wrong
// print some debug information then exit
static inline void _corrupt_paths(const GraphWalker *wlk, size_t num_next,
                                  const dBNode nodes[4],
                                  const Nucleotide bases[4])
{
  size_t i;
  char kstr[MAX_KMER_SIZE+3];

  message("  Fork:\n");

  for(i = 0; i < num_next; i++) {
    db_node_to_str(wlk->db_graph, nodes[i], kstr);
    message("    %s [%c]\n", kstr, dna_nuc_to_char(bases[i]));
  }

  graph_walker_print_state(wlk, ctx_msg_out);

  // Mark next bases available
  bool forks[4] = {0},      taken_curr[4] = {0};
  bool taken_newp[4] = {0}, taken_cntr[4] = {0};

  // mark in d the branches that are available
  for(i = 0; i < num_next; i++) forks[bases[i]] = true;

  // Check for path corruption
  update_path_forks(&wlk->paths,      taken_curr);
  update_path_forks(&wlk->cntr_paths, taken_cntr);

  char bases_fork[20], bases_curr[20], bases_newp[20], bases_cntr[20];
  dna_bases_list_to_str(forks,      bases_fork);
  dna_bases_list_to_str(taken_curr, bases_curr);
  dna_bases_list_to_str(taken_newp, bases_newp);
  dna_bases_list_to_str(taken_cntr, bases_cntr);
  message("forks: %s curr: %s newp: %s cntrp: %s\n",
          bases_fork, bases_curr, bases_newp, bases_cntr);

  warn("Did you build this .ctp against THIS EXACT .ctx? (REALLY?)");
  warn("If you did please report a bug to turner.isaac@gmail.com");
  abort();
}

#define _gw_choose_return(i,s,gap) do {                                     \
  GraphStep _stp = {.idx = (int8_t)(i), .status = (s), .path_gap = (gap)};  \
  return _stp;                                                              \
} while(0)

/**
 * Make a choice at a junction
 * @return index of choice or -1
 */
GraphStep graph_walker_choose(const GraphWalker *wlk, size_t num_next,
                              const dBNode next_nodes[4],
                              const Nucleotide next_bases[4])
{
  #ifdef DEBUG_WALKER
    char kmer_str[MAX_KMER_SIZE+1];
    binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, kmer_str);
    status("  graph_walker_choose(): %s:%i num_next:%zu",
           kmer_str, wlk->node.orient, num_next);
    graph_walker_print_state(wlk, stderr);
  #endif

  const dBGraph *db_graph = wlk->db_graph;

  if(num_next == 0) {
    ctx_assert(wlk->paths.len == 0);
    ctx_assert(wlk->cntr_paths.len == 0);
    _gw_choose_return(-1, GRPHWLK_NOCOVG, 0);
  }

  if(num_next == 1) {
    bool incol = (db_graph->node_in_cols == NULL ||
                  db_node_has_col(db_graph, next_nodes[0].key, wlk->ctxcol));
    _gw_choose_return(0, incol ? GRPHWLK_COLFWD : GRPHWLK_POPFWD, 0);
  }

  int8_t indices[4] = {0,1,2,3};
  dBNode nodes_store[4];
  Nucleotide bases_store[4];
  const dBNode *nodes = next_nodes;
  const Nucleotide* bases = next_bases;
  size_t i, j;

  // Reduce next nodes that are in this colour
  if(db_graph->node_in_cols != NULL)
  {
    nodes = nodes_store;
    bases = bases_store;

    for(i = 0, j = 0; i < num_next; i++)
    {
      if(db_node_has_col(db_graph, next_nodes[i].key, wlk->ctxcol)) {
        nodes_store[j] = next_nodes[i];
        bases_store[j] = next_bases[i];
        indices[j] = (int8_t)i;
        j++;
      }
    }

    num_next = j;

    if(num_next == 0) {
      ctx_assert(wlk->paths.len == 0);
      ctx_assert(wlk->cntr_paths.len == 0);
    }

    if(num_next == 1) _gw_choose_return(indices[0], GRPHWLK_POPFRK_COLFWD, 0);
    if(num_next == 0) _gw_choose_return(-1,         GRPHWLK_NOCOLCOVG,     0);
  }

  // We have hit a fork
  // abandon if no path info
  if(wlk->paths.len == 0) _gw_choose_return(-1, GRPHWLK_NOLINKS, 0);

  // Mark next bases available
  bool forks[4] = {false}, taken[4] = {false};

  // mark in forks the branches that are available
  for(i = 0; i < num_next; i++) forks[bases[i]] = true;

  // Check for path corruption
  update_path_forks(&wlk->paths,      taken);
  update_path_forks(&wlk->cntr_paths, taken);

  if((taken[0] && !forks[0]) || (taken[1] && !forks[1]) ||
     (taken[2] && !forks[2]) || (taken[3] && !forks[3]))
  {
    _corrupt_paths(wlk, num_next, nodes, bases);
  }

  // Do all the oldest paths pick a consistent next node?
  GPathFollow *path, *oldest_path = &wlk->paths.b[0];
  size_t greatest_age;
  Nucleotide greatest_nuc;

  greatest_age = oldest_path->age;
  greatest_nuc = gpath_follow_get_base(oldest_path, oldest_path->pos);

  ctx_assert(oldest_path->pos < oldest_path->len);

  if(greatest_age == 0) _gw_choose_return(-1, GRPHWLK_NOLINKS, 0);

  // Set i to the index of the oldest path to disagree with our oldest path
  // OR wlk->paths.length if all paths agree
  for(i = 1; i < wlk->paths.len; i++) {
    path = &wlk->paths.b[i];
    if(gpath_follow_get_base(path, path->pos) != greatest_nuc) break;
  }

  // If a path of the same age disagrees, cannot proceed
  if(i < wlk->paths.len && wlk->paths.b[i].age == greatest_age)
    _gw_choose_return(-1, GRPHWLK_SPLIT_LINKS, 0);

  size_t choice_age = (i < wlk->paths.len ? wlk->paths.b[i].age : 0);
  const GraphSegment *gseg, *first_seg = mdc_list_getptr(&wlk->gsegs, 0);

  // for(i = 0; i < gseg_list_len(&wlk->gsegs); i++)
  //   printf(" %u", gseg_list_getconstptr(&wlk->gsegs, i)->num_nodes);
  // printf(" [%zu]\n", gseg_list_len(&wlk->gsegs));

  const GraphSegment *choice_seg = first_seg + choice_age;
  while(!choice_seg->in_fork) choice_seg++;

  ctx_assert2(first_seg[greatest_age-1].in_fork, "%zu %u/%u",
              greatest_age, oldest_path->pos, oldest_path->len);
  ctx_assert2(choice_seg < first_seg + greatest_age, "%zu >= %zu",
              choice_seg - first_seg, greatest_age);

  // printf(" choice_seg: %zu\n", choice_seg - first_seg);

  // path_gap is the distance between deciding junctions
  // when this is large we have low confidence
  // see graph_step.h:GraphStep{path_gap}
  size_t path_gap = 0;
  for(gseg = first_seg; gseg <= choice_seg; gseg++)
    path_gap += gseg->num_nodes;

  // Does every next node have a path?
  // Fail if missing assembly info
  if(wlk->missing_path_check &&
     (size_t)taken[0]+taken[1]+taken[2]+taken[3] < num_next) {
    _gw_choose_return(-1, GRPHWLK_MISSING_LINKS, path_gap);
  }

  // There is unique next node
  // Find the correct next node chosen by the paths
  // Assme next node has colour, since we used paths to find it
  //  (paths are colour specify)
  for(i = 0; i < num_next; i++)
    if(bases[i] == greatest_nuc)
      _gw_choose_return(indices[i], GRPHWLK_USELINKS, path_gap);

  // Should be impossible to reach here...
  die("Should be impossible to reach here");
}

#undef _gw_choose_return

/**
 * This is the main traversal function, all other traversal functions call this
 * @param num_nodes is how many nodes we are jumping. If new node is adjacent to
 *                  current node (wlk->node), then num_nodes should be 1.
 * @param lost_nuc  Base lost when moving forward. -1 if moving more than one node
 */
static void _graph_walker_force_jump(GraphWalker *wlk,
                                     dBNode node,
                                     bool is_fork,
                                     size_t num_nodes,
                                     int lost_nuc) // -1 if jump
{
  ctx_assert(node.key != HASH_NOT_FOUND);
  ctx_assert(num_nodes > 0);

  const dBGraph *db_graph = wlk->db_graph;

  #ifdef DEBUG_WALKER
    char kmer_str[MAX_KMER_SIZE+1];
    binary_kmer_to_str(wlk->bkey, db_graph->kmer_size, kmer_str);
    status("  _graph_walker_force_jump(): %s:%i is_fork:%s",
           kmer_str, wlk->node.orient, is_fork ? "yes" : "no");
  #endif

  // last_step should be set by now
  // if we used a path to move to the next node, we must be at a fork
  // The reverse is NOT true, as the caller may have forced us down a junction
  // with graph_walker_force() rather than use a path
  ctx_assert2((wlk->last_step.status != GRPHWLK_USELINKS) || is_fork, "%i vs %i",
              (int)wlk->last_step.status, (int)is_fork);

  // Update GraphWalker position
  wlk->bkey = db_node_get_bkey(db_graph, node.key);
  wlk->node = node;

  if(is_fork)
  {
    // We passed a fork - take all paths that agree with said nucleotide and
    // haven't ended, also update junction progress for each path
    Nucleotide base, pnuc;
    GPathFollow *path;
    size_t i, j, npaths = wlk->paths.len;

    base = bkmer_get_last_nuc(wlk->bkey, node.orient, db_graph->kmer_size);

    // Check curr pathh
    for(i = 0, j = 0; i < npaths; i++)
    {
      path = &wlk->paths.b[i];
      pnuc = gpath_follow_get_base(path, path->pos);
      if(base == pnuc) {
        path->pos++;
        if(path->pos < path->len) {
          wlk->paths.b[j++] = *path;
        }
        else {
          // Finished following a path from start to end
          if(wlk->used_paths) {
            // mark as used
            size_t pathid = gpset_get_pkey(&wlk->gpstore->gpset, path->gpath);
            (void)bitset_set_mt(wlk->used_paths, pathid);
          }
        }
      }
    }

    wlk->paths.len = j;

    // Counter paths
    for(i = 0, j = 0; i < wlk->cntr_paths.len; i++)
    {
      path = &wlk->cntr_paths.b[i];
      pnuc = gpath_follow_get_base(path, path->pos);
      if(base == pnuc && path->pos+1 < path->len) {
        path->pos++;
        wlk->cntr_paths.b[j++] = *path;
      }
    }

    wlk->cntr_paths.len = j;

    // Statistics
    wlk->fork_count++;
  }

  // Find previous nodes
  dBNode prev_nodes[4];
  Nucleotide prev_bases[4];
  size_t num_other_prev = 0;

  if(lost_nuc >= 0 && db_node_in_col(db_graph, wlk->node.key, wlk->ctxcol))
  {
    num_other_prev = db_graph_prev_nodes_with_mask(db_graph, wlk->node,
                                                   (Nucleotide)lost_nuc,
                                                   wlk->ctxcol,
                                                   prev_nodes, prev_bases);

    // Pick up counter paths
    graph_walker_add_counter_paths(wlk, prev_nodes, num_other_prev);
  }

  ctx_assert(!is_fork            || num_nodes == 1);
  ctx_assert(num_other_prev == 0 || num_nodes == 1);

  // Update graph sections
  // num_other_prev > 0 because it doesn't count the node we came from
  _gw_gseg_update(wlk, is_fork, num_other_prev > 0, num_nodes);

  // Pick up new paths
  pickup_paths(wlk, wlk->node, false, 0);
}

/**
 * Jump to a new node within the current sample unitig
 * (can actually be any node up until the end of the current unitig)
 * @param num_nodes is number of nodes we have moved forward
 */
void graph_walker_jump_along_unitig(GraphWalker *wlk, dBNode node, size_t num_nodes)
{
  ctx_assert(num_nodes > 0);

  #ifdef CTXCHECKS
    // This is just a sanity test
    Edges edges = db_node_edges_in_col(db_node_reverse(node), wlk->ctxcol, wlk->db_graph);
    ctx_assert(edges_get_indegree(edges, node.orient) <= 1);
  #endif

  // Need to check if node is in colour
  bool incol = (wlk->db_graph->node_in_cols == NULL ||
                db_node_has_col(wlk->db_graph, node.key, wlk->ctxcol));

  int status = incol ? GRPHWLK_COLFWD : GRPHWLK_POPFWD;
  wlk->last_step = (GraphStep){.idx = 0, .status = status, .path_gap = 0};

  // Don't need to pick up counter paths since there should be none
  // Now do the work
  _graph_walker_force_jump(wlk, node, false, num_nodes, -1);
}

/**
 * Move to the next node
 * @param is_fork If true, node is the result of taking a fork (updates paths)
 */
void graph_walker_force(GraphWalker *wlk, dBNode node, bool is_fork)
{
  ctx_assert(node.key != HASH_NOT_FOUND);
  const size_t kmer_size = wlk->db_graph->kmer_size;
  Nucleotide lost_nuc = bkmer_get_first_nuc(wlk->bkey, wlk->node.orient, kmer_size);

  // _graph_walker_force_jump now picks up counter paths
  _graph_walker_force_jump(wlk, node, is_fork, 1, (int)lost_nuc);
}

// return 1 on success, 0 otherwise
bool graph_walker_next_nodes(GraphWalker *wlk, size_t num_next,
                             const dBNode nodes[4], const Nucleotide bases[4])
{
  wlk->last_step = graph_walker_choose(wlk, num_next, nodes, bases);
  int idx = wlk->last_step.idx;
  if(idx == -1) return false;
  graph_walker_force(wlk, nodes[idx],
                     graph_step_status_is_fork(wlk->last_step.status));
  return true;
}

// return 1 on success, 0 otherwise
bool graph_walker_next(GraphWalker *wlk)
{
  const dBGraph *db_graph = wlk->db_graph;
  Edges edges = db_node_get_edges(db_graph, wlk->node.key, 0);

  dBNode nodes[4];
  Nucleotide bases[4];
  size_t num_next;

  num_next = db_graph_next_nodes(db_graph, wlk->bkey, wlk->node.orient, edges,
                                 nodes, bases);

  return graph_walker_next_nodes(wlk, num_next, nodes, bases);
}


//
// Force traversal along an array of nodes (`priming' a GraphWalker)
//

/**
 * Traversal of every node in a list of nodes using the supplied GraphWalker
 * Visits each node specifed
 */
void graph_walker_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                           bool forward)
{
  Edges edges;
  bool is_fork;
  dBNode next;
  size_t i;
  const dBGraph *db_graph = wlk->db_graph;

  for(i = 0; i < n; i++) {
    edges = db_node_edges_in_col(wlk->node, wlk->ctxcol, db_graph);
    is_fork = edges_get_outdegree(edges, wlk->node.orient) > 1;
    next = db_nodes_get(arr, n, forward, i);
    graph_walker_force(wlk, next, is_fork);
  }
}

void graph_walker_prime(GraphWalker *wlk,
                        const dBNode *block, size_t n,
                        size_t max_context, bool forward)
{
  ctx_assert(n > 0);
  // ctx_check(db_node_check_nodes(block, n, wlk->db_graph));

  // If picking up paths is turned off, jump to last
  if(!gpath_store_use_traverse(wlk->gpstore)) {
    graph_walker_start(wlk, db_nodes_get(block, n, forward, n-1));
    return;
  }

  dBNode node0;

  if(n > max_context) {
    if(forward) block = block + n - max_context;
    n = max_context;
  }

  // DEBUG
  // fprintf(stderr, "\n\n");
  // db_nodes_print(block, n, db_graph, stderr);
  // fprintf(stderr, " orient:%s\n", forward ? "forward" : "reverse");
  //

  if(forward) { node0 = block[0]; block++; }
  else { node0 = db_node_reverse(block[n-1]); }

  graph_walker_start(wlk, node0);
  graph_walker_traverse(wlk, block, n-1, forward);

  // For debugging
  // graph_walker_print_state(wlk, stderr);
}

/**
 * Check the graph walker doesn't veer away from the given contig
 * At each node:
 *  a. If we can't progress -> success
 *  b. If we can and it doesn't match what we expected -> disagrees
 * @param forward Traverse contig forward, otherwise reverse complement nodes
 *                 and work backwards
 */
bool graph_walker_agrees_contig(GraphWalker *wlk,
                                const dBNode *block, size_t num_nodes,
                                bool forward)
{
  if(num_nodes == 0 || !wlk->paths.len) return true;

  size_t i, j, n, njuncs = graph_walker_get_max_path_junctions(wlk);
  dBNode expnode;

  #ifdef CTXCHECKS
    // Check last k-1 bp of bkmer match block
    expnode = db_nodes_get(block, num_nodes, forward, 0);
    BinaryKmer bkmer0 = db_node_oriented_bkmer(wlk->db_graph, wlk->node);
    BinaryKmer bkmer1 = db_node_oriented_bkmer(wlk->db_graph, expnode);
    bkmer0 = binary_kmer_left_shift_one_base(bkmer0, wlk->db_graph->kmer_size);
    binary_kmer_set_last_nuc(&bkmer1, 0);

    if(!binary_kmer_eq(bkmer0, bkmer1))
    {
      char bstr0[MAX_KMER_SIZE+1], bstr1[MAX_KMER_SIZE+1];
      binary_kmer_to_str(bkmer0, wlk->db_graph->kmer_size, bstr0);
      binary_kmer_to_str(bkmer1, wlk->db_graph->kmer_size, bstr1);
      graph_walker_print_state(wlk, stdout);
      printf("wlk: %s contig: %s num_nodes %zu\n", bstr0, bstr1, num_nodes);
    }

    ctx_check(binary_kmer_eq(bkmer0, bkmer1));
  #endif

  dBNode nodes[4];
  Nucleotide nucs[4];
  Edges edges;

  for(i = j = 0; i < num_nodes && j < njuncs; i++, j += (n > 1))
  {
    expnode = db_nodes_get(block, num_nodes, forward, i);
    edges = db_node_get_edges_union(wlk->db_graph, wlk->node.key);

    if(edges_get_outdegree(edges, wlk->node.orient) == 1) {
      nodes[0] = expnode;
      nucs[0] = db_node_get_last_nuc(expnode, wlk->db_graph);
      n = 1;
    } else {
      // Need to look up possible next nodes
      n = db_graph_next_nodes(wlk->db_graph, wlk->bkey, wlk->node.orient,
                              edges, nodes, nucs);
    }

    // If we can't progress -> success
    // if we can and it doesn't match what we expected -> disagrees
    if(!graph_walker_next_nodes(wlk, n, nodes, nucs)) return true;
    if(!db_nodes_are_equal(wlk->node, expnode)) return false;
  }

  return true;
}
