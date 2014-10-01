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

#define USE_COUNTER_PATHS 1

// How many junctions are left to be traversed in our longest remaining path
static size_t graph_walker_get_max_path_junctions(const GraphWalker *wlk)
{
  size_t i, rem, max = 0;
  for(i = 0; i < wlk->paths.len; i++) {
    rem = (size_t)(wlk->paths.data[i].len - wlk->paths.data[i].pos);
    max = MAX2(rem, max);
  }
  return max;
}

static void print_path_list(const GPathFollowBuffer *pbuf, FILE *fout)
{
  size_t i, j;
  GPathFollow *path;

  for(i = 0; i < pbuf->len; i++) {
    path = &pbuf->data[i];
    fprintf(fout, "   %p ", path->gpath->seq);
    for(j = 0; j < path->len; j++)
      fputc(dna_nuc_to_char(gpath_follow_get_base(path, j)), fout);
    fprintf(fout, " [%zu/%zu] age: %zu\n", (size_t)path->pos,
            (size_t)path->len, (size_t)path->age);
  }
}

void graph_walker_print_state(const GraphWalker *wlk, FILE *fout)
{
  char bkmerstr[MAX_KMER_SIZE+1], bkeystr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, bkmerstr);
  binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, bkeystr);
  fprintf(fout, " GWState:%s (%s:%i) ctx: %zu ctp: %zu\n",
         bkmerstr, bkeystr, wlk->node.orient, wlk->ctxcol, wlk->ctpcol);
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

// Only set up memory
// need to call graph_walker_init to reset/initialise state
void graph_walker_alloc(GraphWalker *wlk)
{
  gpath_follow_buf_alloc(&wlk->paths, 256);
  gpath_follow_buf_alloc(&wlk->cntr_paths, 512);
  gseg_list_alloc(&wlk->gsegs, 128);
}

// Free memory
void graph_walker_dealloc(GraphWalker *wlk)
{
  gpath_follow_buf_dealloc(&wlk->paths);
  gpath_follow_buf_dealloc(&wlk->cntr_paths);
  gseg_list_dealloc(&wlk->gsegs);
}

static inline void _gw_gseg_init(GraphWalker *wlk)
{
  ctx_assert(gseg_list_length(&wlk->gsegs) == 0);
  GraphSegment gseg = {.in_fork = 0, .out_fork = 0, .num_nodes = 1};
  gseg_list_unshift(&wlk->gsegs, gseg);
}

static inline void _gw_gseg_update(GraphWalker *wlk,
                                   bool fw_fork, bool rv_fork,
                                   size_t num_nodes)
{
  // First GraphSection is the one the most recent one we saw
  size_t i, len = gseg_list_length(&wlk->gsegs);
  GraphSegment *first = gseg_list_get(&wlk->gsegs, 0);
  ctx_assert(len > 0);

  // Previous section may have ended in a fork - update it
  first->out_fork |= fw_fork;

  if(fw_fork || rv_fork)
  {
    // Start a new section
    ctx_assert(num_nodes == 1);
    GraphSegment gseg = {.in_fork = rv_fork, .out_fork = 0, .num_nodes = 0};
    gseg_list_unshift(&wlk->gsegs, gseg);

    len = gseg_list_length(&wlk->gsegs);
    first = gseg_list_get(&wlk->gsegs, 0);

    // Update all path ages: age is the graph section index
    for(i = 0; i < wlk->paths.len; i++)
      wlk->paths.data[i].age++;

    for(i = 0; i < wlk->cntr_paths.len; i++)
      wlk->cntr_paths.data[i].age++;

    // Drop graph sections without any paths (apart from most recent one)
    size_t max_segs = 1;

    if(wlk->paths.len)
      max_segs = MAX2(max_segs, wlk->paths.data[0].age + 1);
    if(wlk->cntr_paths.len)
      max_segs = MAX2(max_segs, wlk->cntr_paths.data[0].age + 1);

    ctx_assert2(max_segs <= len, "%zu > %zu", max_segs, len);
    gseg_list_popn(&wlk->gsegs, len - max_segs);
  }

  // printf("jump %zu\n", num_nodes);
  first->num_nodes += num_nodes;
}

// Returns number of paths picked up
// next_nuc only used if counter == true and node has out-degree > 1
static inline size_t pickup_paths(GraphWalker *wlk, dBNode node,
                                  bool counter, Nucleotide next_nuc)
{
  // printf("pickup %s paths from: %zu:%i\n", counter ? "counter" : "curr",
  //        (size_t)index, orient);

  const dBGraph *db_graph = wlk->db_graph;
  const GPathStore *gpstore = wlk->gpstore;
  GPathFollowBuffer *pbuf = counter ? &wlk->cntr_paths : &wlk->paths;
  const size_t ncols = gpstore->gpset.ncols, num_paths = pbuf->len;

  // Picking up paths is turned off
  if(!gpath_store_use_traverse(gpstore)) return 0;

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
        gpath_follow_buf_add(pbuf, fpath);
      }
    }
  }

  #ifdef DEBUG_WALKER
    char bkey_str[MAX_KMER_SIZE+1], node_str[MAX_KMER_SIZE+1];
    BinaryKmer node_bkey = db_node_get_bkmer(wlk->db_graph, node.key);
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
 * @param prev_nodes nodes before wlk->node, oriented away from wlk->node
 *        i.e. the nodes you would reach if you walked from
 *        db_node_reverse(wlk->node)
*/
void graph_walker_add_counter_paths(GraphWalker *wlk,
                                    const dBNode prev_nodes[4],
                                    size_t num_prev)
{
  size_t i;
  Nucleotide next_base = binary_kmer_last_nuc(wlk->bkmer);

  // Reverse orientation, pick up paths
  for(i = 0; i < num_prev; i++)
    pickup_paths(wlk, db_node_reverse(prev_nodes[i]), true, next_base);
}


void graph_walker_init(GraphWalker *wlk, const dBGraph *graph,
                       Colour ctxcol, Colour ctpcol, dBNode node)
{
  // Check that the graph is loaded properly (all edges merged into one colour)
  ctx_assert(graph->num_edge_cols == 1);
  ctx_assert(graph->num_of_cols == 1 || graph->node_in_cols != NULL);

  ctx_assert(wlk->paths.len == 0);
  ctx_assert(wlk->cntr_paths.len == 0);
  ctx_assert(gseg_list_length(&wlk->gsegs) == 0);

  wlk->db_graph = graph;
  wlk->gpstore = &graph->gpstore;
  wlk->ctxcol = ctxcol;
  wlk->ctpcol = ctpcol,
  wlk->node = node;

  // stats
  wlk->fork_count = 0;
  wlk->last_step = (GraphStep){.idx = -1, 0};

  // Get bkmer oriented correctly (not bkey)
  wlk->bkey = db_node_get_bkmer(graph, node.key);
  wlk->bkmer = db_node_oriented_bkmer(graph, node);

  #ifdef DEBUG_WALKER
    char kmer_str[MAX_KMER_SIZE+1];
    binary_kmer_to_str(wlk->bkey, graph->kmer_size, kmer_str);
    status("  graph_walker_init(): %s:%i", kmer_str, wlk->node.orient);
  #endif

  _gw_gseg_init(wlk);

  // Pick up new paths
  pickup_paths(wlk, wlk->node, false, 0);

  ctx_assert(gseg_list_length(&wlk->gsegs) == 1);
  // printf("zero: %u\n", gseg_list_get(&wlk->gsegs, 0)->num_nodes);

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
  hash = CityHash64WithSeeds((const char*)wlk->paths.data, path_bytes,
                             hash, wlk->paths.len);
  hash = CityHash64WithSeeds((const char*)wlk->cntr_paths.data, cntr_path_bytes,
                             hash, wlk->cntr_paths.len);

  return hash;
}

//
// Junction decision
//

static inline void update_path_forks(GPathFollowBuffer *pbuf, bool taken[4])
{
  size_t i;
  GPathFollow *path;
  for(i = 0; i < pbuf->len; i++) {
    path = &pbuf->data[i];
    taken[gpath_follow_get_base(path, path->pos)] = true;
  }
}

static inline void _corrupt_paths(GraphWalker *wlk, size_t num_next,
                                  const dBNode nodes[4],
                                  const Nucleotide bases[4])
__attribute__((noreturn));

// If we call this funcitons, something has gone wrong
// print some debug information then exit
static inline void _corrupt_paths(GraphWalker *wlk, size_t num_next,
                                  const dBNode nodes[4],
                                  const Nucleotide bases[4])
{
  size_t i;
  BinaryKmer bkey;
  char str[MAX_KMER_SIZE+1];

  message("  Fork:\n");

  for(i = 0; i < num_next; i++) {
    bkey = db_node_get_bkmer(wlk->db_graph, nodes[i].key);
    binary_kmer_to_str(bkey, wlk->db_graph->kmer_size, str);
    message("    %s:? [%c]\n", str, dna_nuc_to_char(bases[i]));
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

#define _gw_choose_return(i,s,hascol,gap) do { \
  GraphStep _stp = {.idx = (int8_t)(i),        \
                    .status = (s),             \
                    .node_has_col = (hascol),  \
                    .path_gap = (gap)          \
                   };                          \
  return _stp;                                 \
} while(0)

/**
 * Make a choice at a junction
 * GraphWalker is not const because we update the path junction cache
 * @return index of choice or -1
 */
GraphStep graph_walker_choose(GraphWalker *wlk, size_t num_next,
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

  if(num_next == 0) {
    if(wlk->paths.len > 0) {
      graph_walker_print_state(wlk, stderr);
      die("Shouldn't be able to reach here");
    }

    ctx_assert(wlk->paths.len == 0);
    ctx_assert(wlk->cntr_paths.len == 0);
    _gw_choose_return(-1, GRPHWLK_NOCOVG, false, 0);
  }

  const dBGraph *db_graph = wlk->db_graph;

  if(num_next == 1) {
    bool incol = (db_graph->node_in_cols == NULL ||
                  db_node_has_col(db_graph, next_nodes[0].key, wlk->ctxcol));
    _gw_choose_return(0, GRPHWLK_FORWARD, incol, 0);
  }

  int8_t indices[4] = {0,1,2,3};
  dBNode nodes_store[4];
  Nucleotide bases_store[4];
  const dBNode *nodes = nodes_store;
  const Nucleotide* bases = bases_store;
  size_t i, j;

  // Reduce next nodes that are in this colour
  if(db_graph->node_in_cols != NULL)
  {
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

    if(num_next == 1) _gw_choose_return(indices[0], GRPHWLK_COLFWD,    true, 0);
    if(num_next == 0) _gw_choose_return(-1,         GRPHWLK_NOCOLCOVG, false,0);
  }
  else {
    nodes = next_nodes;
    bases = next_bases;
  }

  // We have hit a fork
  // abandon if no path info
  if(wlk->paths.len == 0) _gw_choose_return(-1, GRPHWLK_NOPATHS, false, 0);

  // Mark next bases available
  bool forks[4] = {false}, taken[4] = {false};

  // mark in d the branches that are available
  for(i = 0; i < num_next; i++) forks[bases[i]] = 1;

  // Check for path corruption
  update_path_forks(&wlk->paths,      taken);
  update_path_forks(&wlk->cntr_paths, taken);

  if((taken[0] && !forks[0]) || (taken[1] && !forks[1]) ||
     (taken[2] && !forks[2]) || (taken[3] && !forks[3])) {
    _corrupt_paths(wlk, num_next, nodes, bases);
  }

  // Do all the oldest paths pick a consistent next node?
  GPathFollow *path, *oldest_path = &wlk->paths.data[0];
  size_t greatest_age;
  Nucleotide greatest_nuc;

  greatest_age = oldest_path->age;
  greatest_nuc = gpath_follow_get_base(oldest_path, oldest_path->pos);

  ctx_assert(oldest_path->pos < oldest_path->len);

  if(greatest_age == 0) _gw_choose_return(-1, GRPHWLK_NOPATHS, false, 0);

  // Set i to the index of the oldest path to disagree with our oldest path
  // OR wlk->paths.length if all paths agree
  for(i = 1; i < wlk->paths.len; i++) {
    path = &wlk->paths.data[i];
    if(gpath_follow_get_base(path, path->pos) != greatest_nuc) break;
  }

  // If a path of the same age disagrees, cannot proceed
  if(i < wlk->paths.len && wlk->paths.data[i].age == greatest_age)
    _gw_choose_return(-1, GRPHWLK_SPLIT_PATHS, false, 0);

  size_t choice_age = (i < wlk->paths.len ? wlk->paths.data[i].age : 0);
  GraphSegment *gseg, *choice_seg, *first_seg = gseg_list_get(&wlk->gsegs, 0);

  // for(i = 0; i < gseg_list_length(&wlk->gsegs); i++)
  //   printf(" %u", gseg_list_get(&wlk->gsegs, i)->num_nodes);
  // printf(" [%zu]\n", gseg_list_length(&wlk->gsegs));

  choice_seg = first_seg + choice_age;
  while(!choice_seg->in_fork) choice_seg++;

  ctx_assert2(first_seg[greatest_age-1].in_fork, "%zu %u/%u", greatest_age, oldest_path->pos, oldest_path->len);
  ctx_assert2(choice_seg < first_seg + greatest_age, "%zu >= %zu",
              choice_seg - first_seg, greatest_age);

  // printf(" choice_seg: %zu\n", choice_seg - first_seg);

  // path_gap is the distance between deciding junctions
  // when this is large we have low confidence
  // see graph_step.h:GraphStep{path_gap}
  size_t path_gap = 0;
  for(gseg = first_seg; gseg <= choice_seg; gseg++)
    path_gap += gseg->num_nodes;

  #ifdef USE_COUNTER_PATHS
  // Does every next node have a path?
  // Fail if missing assembly info
  if((size_t)taken[0]+taken[1]+taken[2]+taken[3] < num_next)
    _gw_choose_return(-1, GRPHWLK_MISSING_PATHS, false, path_gap);
  #endif

  // There is unique next node
  // Find the correct next node chosen by the paths
  // Assme next node has colour, since we used paths to find it
  //  (paths are colour specify)
  for(i = 0; i < num_next; i++)
    if(bases[i] == greatest_nuc)
      _gw_choose_return(indices[i], GRPHWLK_USEPATH, true, path_gap);

  // Should be impossible to reach here...
  die("Should be impossible to reach here");
}

#undef _gw_choose_return

// This is the main traversal function, all other traversal functions call this
// in_degree_fork indicates if in degree of new node is greater than one
//  (in this colour)
static void _graph_walker_force_jump(GraphWalker *wlk,
                                     hkey_t hkey, BinaryKmer bkmer,
                                     bool is_fork,
                                     size_t num_nodes,
                                     int lost_nuc) // -1 if jump
{
  ctx_assert(hkey != HASH_NOT_FOUND);

  #ifdef DEBUG_WALKER
    char kmer_str[MAX_KMER_SIZE+1];
    binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, kmer_str);
    status("  _graph_walker_force_jump(): %s:%i is_fork:%s",
           kmer_str, wlk->node.orient, is_fork ? "yes" : "no");
  #endif

  if(is_fork)
  {
    // We passed a fork - take all paths that agree with said nucleotide and
    // haven't ended, also update junction progress for each path
    Nucleotide base = binary_kmer_last_nuc(bkmer), pnuc;
    GPathFollow *path;
    size_t i, j, npaths = wlk->paths.len;

    // Check curr pathh
    for(i = 0, j = 0; i < npaths; i++)
    {
      path = &wlk->paths.data[i];
      pnuc = gpath_follow_get_base(path, path->pos);
      if(base == pnuc && path->pos+1 < path->len) {
        path->pos++;
        wlk->paths.data[j++] = *path;
      }
    }

    wlk->paths.len = j;

    // Counter paths
    for(i = 0, j = 0; i < wlk->cntr_paths.len; i++)
    {
      path = &wlk->cntr_paths.data[i];
      pnuc = gpath_follow_get_base(path, path->pos);
      if(base == pnuc && path->pos+1 < path->len) {
        path->pos++;
        wlk->cntr_paths.data[j++] = *path;
      }
    }

    wlk->cntr_paths.len = j;

    // Statistics
    wlk->fork_count++;
  }

  // Update GraphWalker position
  wlk->bkmer = bkmer;
  wlk->bkey = db_node_get_bkmer(wlk->db_graph, hkey);
  wlk->node.key = hkey;
  wlk->node.orient = bkmer_get_orientation(wlk->bkmer, wlk->bkey);

  // Find previous nodes
  dBNode prev_nodes[4];
  Nucleotide prev_bases[4];
  size_t num_prev = 0;

  if(lost_nuc >= 0)
  {
    num_prev = db_graph_prev_nodes_with_mask(wlk->db_graph, wlk->node,
                                             (Nucleotide)lost_nuc,
                                             wlk->ctxcol,
                                             prev_nodes, prev_bases);

    // Pick up counter paths
    graph_walker_add_counter_paths(wlk, prev_nodes, num_prev);
  }

  ctx_assert(!is_fork      || num_nodes == 1);
  ctx_assert(num_prev == 0 || num_nodes == 1);

  // Update graph sections
  _gw_gseg_update(wlk, is_fork, num_prev > 0, num_nodes);

  // Pick up new paths
  pickup_paths(wlk, wlk->node, false, 0);
}

/**
 * Jump to a new node within the current sample supernode
 * (can actually be any node up until the end of the current supernode)
 * @param num_nodes is number of nodes we have moved forward
 */
void graph_walker_jump_along_snode(GraphWalker *wlk, hkey_t hkey,
                                   BinaryKmer bkmer, size_t num_nodes)
{
  #ifdef CTXCHECKS
    // This is just a sanity test
    Edges edges = db_node_get_edges(wlk->db_graph, hkey, 0);
    BinaryKmer bkey = db_node_get_bkmer(wlk->db_graph, hkey);
    Orientation orient = bkmer_get_orientation(bkmer, bkey);
    ctx_assert(edges_get_indegree(edges, orient) <= 1);
  #endif

  // Don't need to pick up counter paths since there should be none
  // Now do the work
  _graph_walker_force_jump(wlk, hkey, bkmer, false, num_nodes, -1);
}

/**
 * Move to the next node
 */
void graph_walker_force(GraphWalker *wlk, hkey_t hkey, Nucleotide base,
                        bool is_fork)
{
  ctx_assert(hkey != HASH_NOT_FOUND);
  BinaryKmer bkmer;
  const size_t kmer_size = wlk->db_graph->kmer_size;
  Nucleotide lost_nuc = binary_kmer_first_nuc(wlk->bkmer, kmer_size);
  bkmer = binary_kmer_left_shift_add(wlk->bkmer, kmer_size, base);

  // _graph_walker_force_jump now picks up counter paths
  _graph_walker_force_jump(wlk, hkey, bkmer, is_fork, 1, (int)lost_nuc);
}

// return 1 on success, 0 otherwise
bool graph_walker_next_nodes(GraphWalker *wlk, size_t num_next,
                             const dBNode nodes[4], const Nucleotide bases[4])
{
  wlk->last_step = graph_walker_choose(wlk, num_next, nodes, bases);
  int idx = wlk->last_step.idx;
  if(idx == -1) return false;
  graph_walker_force(wlk, nodes[idx].key, bases[idx],
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

// Fast traverse - avoid a bkmer_revcmp
static inline void _graph_walker_fast(GraphWalker *wlk,
                                      const dBNode next_node,
                                      bool is_fork,
                                      size_t num_nodes)
{
  const size_t kmer_size = wlk->db_graph->kmer_size;
  BinaryKmer bkmer, bkey;
  Nucleotide nuc;
  ctx_assert(num_nodes > 0);

  // Only one path between two nodes
  if(num_nodes == 1) {
    nuc = db_node_get_last_nuc(next_node, wlk->db_graph);
    graph_walker_force(wlk, next_node.key, nuc, is_fork);
  }
  else {
    // jumping to the end of a supernode
    bkey = db_node_get_bkmer(wlk->db_graph, next_node.key);
    bkmer = bkmer_oriented_bkmer(bkey, next_node.orient, kmer_size);
    graph_walker_jump_along_snode(wlk, next_node.key, bkmer, num_nodes);
  }

  // char tmpbkmer[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, tmpbkmer);
  // printf("  forced: %s\n", tmpbkmer);
}

// DEV: not currently working
// Fast traversal of a list of nodes using the supplied GraphWalker
// Only visits nodes deemed informative + last node
// Must have previously initialised or walked to the prior node,
// using: graph_walker_init, graph_walker_force, graph_walker_jump_along_snode,
// graph_traverse or graph_walker_next_nodes
// i.e. wlk->node is a node adjacent to arr[0]
void graph_walker_fast_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                bool forward)
{
  if(n == 0) return;

  // Only one colour should be loaded
  // (so we don't have to use union edges to figure out forks in given colour)
  ctx_assert(wlk->db_graph->num_edge_cols == 1);

  size_t i, num_nodes, walk_start = 0;
  bool infork[3] = {false, false, false}, outfork[3] = {false, false, false};
  Edges edges;
  dBNode nodes[3];

  edges = db_node_get_edges(wlk->db_graph, wlk->node.key, 0);
  outfork[0] = edges_get_outdegree(edges, wlk->node.orient) > 1;

  nodes[0] = wlk->node;
  nodes[1] = forward ? arr[0] : db_node_reverse(arr[n-1]);

  edges = db_node_get_edges(wlk->db_graph, nodes[1].key, 0);
  outfork[1] = edges_get_outdegree(edges, nodes[1].orient) > 1;
  infork[1] = edges_get_indegree(edges, nodes[1].orient) > 1;

  for(i = 0; i+1 < n; i++)
  {
    // Move to node i if informative (given nodes i-1, i+1)
    // node i refers to infork[1] and outfork[1]

    nodes[2] = forward ? arr[i+1] : db_node_reverse(arr[n-i-2]);

    edges = db_node_get_edges(wlk->db_graph, nodes[2].key, 0);
    outfork[2] = edges_get_outdegree(edges, nodes[2].orient) > 1;
    infork[2] = edges_get_indegree(edges, nodes[2].orient) > 1;

    // Traverse nodes[i] if:
    // - previous node had out-degree > 1 (update/drop paths)
    // - current node has in-degree > 1 (pick up counter-paths + merge in new paths)
    // - next node has in-degree > 1 (pickup paths)
    if(outfork[0] || infork[1] || infork[2]) {
      num_nodes = i - walk_start + 1;
      _graph_walker_fast(wlk, nodes[1], outfork[0], num_nodes);
      walk_start = i + 1;
    }

    // Rotate edges, nodes
    infork[0] = infork[1]; infork[1] = infork[2];
    outfork[0] = outfork[1]; outfork[1] = outfork[2];
    nodes[0] = nodes[1]; nodes[1] = nodes[2];
  }

  // Traverse last node
  num_nodes = n - walk_start;
  _graph_walker_fast(wlk, nodes[1], outfork[0], num_nodes);
}

// Traversal of every node in a list of nodes using the supplied GraphWalker
// Visits each node specifed
void graph_walker_slow_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                bool forward)
{
  Edges edges;
  bool is_fork;
  dBNode next;
  Nucleotide nuc;
  size_t i;
  const dBGraph *db_graph = wlk->db_graph;

  for(i = 0; i < n; i++) {
    edges = db_node_get_edges(db_graph, wlk->node.key, 0);
    is_fork = edges_get_outdegree(edges, wlk->node.orient) > 1;
    next = forward ? arr[i] : db_node_reverse(arr[n-1-i]);
    nuc = db_node_get_last_nuc(next, db_graph);
    graph_walker_force(wlk, next.key, nuc, is_fork);
  }
}

void graph_walker_prime(GraphWalker *wlk,
                        const dBNode *block, size_t n,
                        size_t max_context, bool forward,
                        size_t ctxcol, size_t ctpcol,
                        const dBGraph *db_graph)
{
  ctx_assert(n > 0);
  ctx_check(db_node_check_nodes(block, n, db_graph));

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

  graph_walker_init(wlk, db_graph, ctxcol, ctpcol, node0);
  // graph_walker_fast_traverse(wlk, block, n-1, forward);
  graph_walker_slow_traverse(wlk, block, n-1, forward);

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
  dBNode nodes[4], expnode;
  Nucleotide nucs[4];
  Edges edges;

  #ifdef CTXCHECKS
    // Check last k-1 bp of wlk->bkmer match block
    expnode = forward ? block[0] : db_node_reverse(block[num_nodes-1]);
    BinaryKmer bkmer0 = wlk->bkmer;
    BinaryKmer bkmer1 = db_node_oriented_bkmer(wlk->db_graph, expnode);
    bkmer0 = binary_kmer_left_shift_one_base(bkmer0, wlk->db_graph->kmer_size);
    binary_kmer_set_last_nuc(&bkmer1, 0);

    if(!binary_kmers_are_equal(bkmer0, bkmer1))
    {
      char bstr0[MAX_KMER_SIZE+1], bstr1[MAX_KMER_SIZE+1];
      binary_kmer_to_str(bkmer0, wlk->db_graph->kmer_size, bstr0);
      binary_kmer_to_str(bkmer1, wlk->db_graph->kmer_size, bstr1);
      graph_walker_print_state(wlk, stdout);
      printf("wlk: %s contig: %s num_nodes %zu\n", bstr0, bstr1, num_nodes);
    }

    ctx_check(binary_kmers_are_equal(bkmer0, bkmer1));
  #endif

  for(i = 0, j = 0; i < num_nodes && j < njuncs; i++, j += (n > 1))
  {
    expnode = forward ? block[i] : db_node_reverse(block[num_nodes-i-1]);
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
