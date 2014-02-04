#include "global.h"
#include "util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_walker.h"
#include "packed_path.h"

// hash functions
#include "jenkins.h"
#include "twang.h"

#ifdef CTXVERBOSE
#define DEBUG_WALKER 1
#endif

#define USE_COUNTER_PATHS 1

// For GraphWalker to work we assume all edges are merged into one colour
// (i.e. graph->num_edge_cols == 1)
// If only one colour loaded we assume all edges belong to this colour

FollowPath follow_path_create(const uint8_t *seq, PathLen plen)
{
  FollowPath fpath = {.seq = seq, .len = plen, .pos = 0, .first_cached = 0};
  memcpy(fpath.cache, fpath.seq, sizeof(fpath.cache));
  return fpath;
}

// Get a base from the FollowPath cache
static inline Nucleotide cache_fetch(const FollowPath *path, size_t pos)
{
  return packed_fetch(path->cache, pos - path->first_cached);
}

// Check if the FollowPath cache needs updated, based of path->pos value
// if it does, update it
static inline void cache_update(FollowPath *path)
{
  // 4 bases per byte
  PathLen new_cache_start = sizeof(path->cache) * 4 *
                            (path->pos / (sizeof(path->cache) * 4));

  if(new_cache_start != path->first_cached) {
    path->first_cached = new_cache_start;
    memcpy(path->cache, path->seq + path->first_cached/4, sizeof(path->cache));
  }
}

static void print_path_list(const PathBuffer *pbuf)
{
  size_t i, j;
  FollowPath *path;

  for(i = 0; i < pbuf->len; i++) {
    path = &pbuf->data[i];
    printf("   %p ", path->seq);
    for(j = 0; j < path->len; j++)
      putc(dna_nuc_to_char(cache_fetch(path, j)), stdout);
    printf(" [%zu/%zu]\n", (size_t)path->pos, (size_t)path->len);
  }
}

void graph_walker_print_state(const GraphWalker *wlk)
{
  char bkmerstr[MAX_KMER_SIZE+1], bkeystr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, bkmerstr);
  binary_kmer_to_str(wlk->bkey, wlk->db_graph->kmer_size, bkeystr);
  printf(" GWState:%s (%s:%i) ctx: %zu ctp: %zu\n",
         bkmerstr, bkeystr, wlk->node.orient, wlk->ctxcol, wlk->ctpcol);
  printf("  num_curr: %zu\n", wlk->paths.len);
  print_path_list(&wlk->paths);
  printf("  num_new: %zu\n", wlk->new_paths.len);
  print_path_list(&wlk->new_paths);
  printf("  num_counter: %zu\n", wlk->cntr_paths.len);
  print_path_list(&wlk->cntr_paths);
  printf("--\n");
}

size_t graph_walker_est_mem()
{
  return sizeof(FollowPath)*1024;
}

// Only set up memory
// need to call graph_walker_init to reset/initialise state
void graph_walker_alloc(GraphWalker *wlk)
{
  path_buf_alloc(&wlk->paths, 256);
  path_buf_alloc(&wlk->new_paths, 128);
  path_buf_alloc(&wlk->cntr_paths, 512);
}

// Free memory
void graph_walker_dealloc(GraphWalker *wlk)
{
  path_buf_dealloc(&wlk->paths);
  path_buf_dealloc(&wlk->new_paths);
  path_buf_dealloc(&wlk->cntr_paths);
}

// Returns number of paths picked up
// next_nuc only used if counter == true and node has out-degree > 1
static inline size_t pickup_paths(GraphWalker *wlk, dBNode node,
                                  boolean counter, Nucleotide next_nuc)
{
  // printf("pickup %s paths from: %zu:%i\n", counter ? "counter" : "curr",
  //        (size_t)index, orient);

  const dBGraph *db_graph = wlk->db_graph;
  const PathStore *pstore = wlk->pstore;
  PathBuffer *pbuf = counter ? &wlk->cntr_paths : &wlk->new_paths;
  size_t num_paths = pbuf->len;
  PathIndex pindex;
  PathLen plen;
  Orientation porient;
  const uint8_t *path, *seq;
  boolean cntr_filter_nuc0 = false;

  if(counter) {
    // cntr_filter_nuc0 is needed for picking up counter paths with outdegree > 1
    cntr_filter_nuc0 = (db_node_outdegree_in_col(node, wlk->ctxcol, db_graph) > 1);
  }

  pindex = *db_node_paths_volptr(db_graph, node.key);

  while(pindex != PATH_NULL)
  {
    path = pstore->store+pindex;
    packedpath_get_len_orient(path, pstore->colset_bytes, &plen, &porient);

    if(node.orient == porient && packedpath_has_col(path, wlk->ctpcol))
    {
      seq = packedpath_seq(path, pstore->colset_bytes);
      FollowPath fpath = follow_path_create(seq, plen);

      if(!cntr_filter_nuc0 || cache_fetch(&fpath, 0) == next_nuc) {
        if(cntr_filter_nuc0) fpath.pos++; // already took a base

        path_buf_add(pbuf, fpath);
      }
    }

    pindex = packedpath_get_prev(pstore->store+pindex);
  }

  return pbuf->len - num_paths;
}

void graph_walker_init(GraphWalker *wlk, const dBGraph *graph,
                       Colour ctxcol, Colour ctpcol, dBNode node)
{
  // Check that the graph is loaded properly (all edges merged into one colour)
  assert(graph->num_edge_cols == 1);
  assert(graph->num_of_cols == 1 || graph->node_in_cols != NULL);

  GraphWalker gw = {.db_graph = graph, .pstore = &graph->pdata,
                    .ctxcol = ctxcol, .ctpcol = ctpcol,
                    .node = node,
                    // new
                    .paths = wlk->paths, .cntr_paths = wlk->cntr_paths,
                    // stats
                    .fork_count = 0, .last_step = {.idx = -1, 0}};

  memcpy(wlk, &gw, sizeof(GraphWalker));

  path_buf_reset(&wlk->paths);
  path_buf_reset(&wlk->new_paths);
  path_buf_reset(&wlk->cntr_paths);

  // Get bkmer oriented correctly (not bkey)
  wlk->bkey = db_node_get_bkmer(graph, node.key);
  wlk->bkmer = db_graph_oriented_bkmer(graph, node.key, node.orient);

  // Pick up new paths
  pickup_paths(wlk, wlk->node, false, 0);
}

void graph_walker_finish(GraphWalker *wlk)
{
  path_buf_reset(&wlk->paths);
  path_buf_reset(&wlk->new_paths);
  path_buf_reset(&wlk->cntr_paths);
}

//
// Hash function
//

// Hash a path using its length, sequence and current offset/position
static inline uint32_t follow_path_hash(const FollowPath *path, uint32_t hash32)
{
  // 16 bases per 32 bit word
  size_t i; uint8_t key[4];
  memcpy(key, &path->len, 2);
  memcpy(key+2, &path->pos, 2);
  hash32 = jenkins_mix(hash32, key[0]);
  hash32 = jenkins_mix(hash32, key[1]);
  hash32 = jenkins_mix(hash32, key[2]);
  hash32 = jenkins_mix(hash32, key[3]);

  for(i = 0; i < sizeof(path->cache); i++)
    hash32 = jenkins_mix(hash32, path->cache[i]);

  // note: we done call jenkins_finish(hash32)
  return hash32;
}

static inline uint32_t follow_pathbuf_hash(const PathBuffer *pbuf, uint32_t hash32)
{
  size_t i;

  for(i = 0; i < pbuf->len; i++)
    hash32 = follow_path_hash(&pbuf->data[i], hash32);

  // note: we done call jenkins_finish(hash32)
  return hash32;
}

// Hash a binary kmer + GraphWalker paths with offsets
uint32_t graph_walker_hash(const GraphWalker *wlk)
{
  // reduce to node from 64 -> 32 bits
  uint64_t hash64 = (wlk->node.key<<1) | wlk->node.orient;
  uint_fast32_t hash32 = twang_hash64to32(hash64);

  // mix in path information
  hash32 = follow_pathbuf_hash(&wlk->paths, hash32);
  hash32 = follow_pathbuf_hash(&wlk->new_paths, hash32);
  hash32 = follow_pathbuf_hash(&wlk->cntr_paths, hash32);

  return jenkins_finish(hash32);
}

//
// Junction decision
//


static inline void update_path_forks(const PathBuffer *pbuf, uint8_t taken[4])
{
  size_t i;
  const FollowPath *path;
  for(i = 0; i < pbuf->len; i++) {
    path = &pbuf->data[i];
    taken[cache_fetch(path, path->pos)] = 1;
  }
}

// out must be at least 11 bytes long: "A, C, G, T"
static inline size_t bases_list_to_str(const uint8_t bases[4], char *out)
{
  size_t i;
  char *str = out;
  const char seq[] = "ACGT";
  for(i = 0; i < 4; i++) {
    if(bases[i]) {
      if(str > out) { memcpy(str, ", ", 2); str += 2; }
      *str = seq[i];
      str++;
    }
  }
  *str = '\0';
  return str-out;
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
  BinaryKmer bkey;
  char str[MAX_KMER_SIZE+1];

  message("  Fork:\n");

  for(i = 0; i < num_next; i++) {
    bkey = db_node_get_bkmer(wlk->db_graph, nodes[i].key);
    binary_kmer_to_str(bkey, wlk->db_graph->kmer_size, str);
    message("    %s:? [%c]\n", str, dna_nuc_to_char(bases[i]));
  }

  graph_walker_print_state(wlk);


  // Mark next bases available
  uint8_t forks[4], taken_curr[4], taken_newp[4], taken_cntr[4];
  memset(forks, 0, sizeof(forks));
  memset(taken_curr, 0, sizeof(taken_curr));
  memset(taken_newp, 0, sizeof(taken_newp));
  memset(taken_cntr, 0, sizeof(taken_cntr));

  // mark in d the branches that are available
  for(i = 0; i < num_next; i++) forks[bases[i]] = 1;

  // Check for path corruption
  update_path_forks(&wlk->paths, taken_curr);
  update_path_forks(&wlk->new_paths, taken_newp);
  update_path_forks(&wlk->cntr_paths, taken_cntr);

  char bases_fork[20], bases_curr[20], bases_newp[20], bases_cntr[20];
  bases_list_to_str(forks, bases_fork);
  bases_list_to_str(taken_curr, bases_curr);
  bases_list_to_str(taken_newp, bases_newp);
  bases_list_to_str(taken_cntr, bases_cntr);
  message("forks: %s curr: %s newp: %s cntrp: %s\n",
          bases_fork, bases_curr, bases_newp, bases_cntr);

  warn("Did you build this .ctp against THIS EXACT .ctx? (REALLY?)");
  abort();
}

#define return_step(i,s) { \
  GraphStep _stp = {.idx = (int8_t)(i), .status = (s)}; \
  return _stp; \
}

// Returns index of choice or -1
// Sets is_fork_in_col true if there is a fork in the given colour
GraphStep graph_walker_choose(const GraphWalker *wlk, size_t num_next,
                              const dBNode next_nodes[4],
                              const Nucleotide next_bases[4])
{
  // #ifdef DEBUG_WALKER
  //   printf("CHOOSE\n");
  //   print_state(wlk);
  // #endif

  if(num_next == 0) {
    assert(wlk->paths.len == 0);
    assert(wlk->new_paths.len == 0);
    assert(wlk->cntr_paths.len == 0);
  }

  if(num_next == 0) return_step(-1, GRPHWLK_NOCOVG);
  if(num_next == 1) return_step( 0, GRPHWLK_FORWARD);

  int8_t indices[4] = {0,1,2,3};
  dBNode nodes_store[4];
  Nucleotide bases_store[4];
  const dBNode *nodes = nodes_store;
  const Nucleotide* bases = bases_store;
  size_t i, j;

  // Reduce next nodes that are in this colour
  if(wlk->db_graph->num_of_cols > 1)
  {
    for(i = 0, j = 0; i < num_next; i++)
    {
      if(db_node_has_col(wlk->db_graph, next_nodes[i].key, wlk->ctxcol)) {
        nodes_store[j] = next_nodes[i];
        bases_store[j] = next_bases[i];
        indices[j] = (int8_t)i;
        j++;
      }
    }

    num_next = j;

    if(num_next == 0) {
      assert(wlk->paths.len == 0);
      assert(wlk->new_paths.len == 0);
      assert(wlk->cntr_paths.len == 0);
    }

    if(num_next == 1) return_step(indices[0], GRPHWLK_COLFWD);
    if(num_next == 0) return_step(-1,         GRPHWLK_NOCOLCOVG);
  }
  else {
    nodes = next_nodes;
    bases = next_bases;
  }

  // We have hit a fork
  // abandon if no path info
  if(wlk->paths.len == 0) return_step(-1, GRPHWLK_NOPATHS);

  // Mark next bases available
  uint8_t forks[4] = {0,0,0,0}, taken[4] = {0,0,0,0};

  // mark in d the branches that are available
  for(i = 0; i < num_next; i++) forks[bases[i]] = 1;

  // Check for path corruption
  update_path_forks(&wlk->paths, taken);
  update_path_forks(&wlk->new_paths, taken);
  update_path_forks(&wlk->cntr_paths, taken);

  if((taken[0] && !forks[0]) || (taken[1] && !forks[1]) ||
     (taken[2] && !forks[2]) || (taken[3] && !forks[3])) {
    _corrupt_paths(wlk, num_next, nodes, bases);
  }

  // Do all the oldest paths pick a consistent next node?
  FollowPath *path, *oldest_path = &wlk->paths.data[0];
  PathLen greatest_age;
  Nucleotide greatest_nuc;

  greatest_age = oldest_path->pos;
  greatest_nuc = cache_fetch(oldest_path, oldest_path->pos);

  for(i = 1; i < wlk->paths.len; i++) {
    path = &wlk->paths.data[i];
    if(path->pos < greatest_age) break;
    if(cache_fetch(path, path->pos) != greatest_nuc)
      return_step(-1, GRPHWLK_SPLIT_PATHS);
  }

  #ifdef USE_COUNTER_PATHS
  // Does every next node have a path?
  // Fail if missing assembly info
  if((size_t)taken[0]+taken[1]+taken[2]+taken[3] < num_next)
    return_step(-1, GRPHWLK_MISSING_PATHS);
  #endif

  // There is unique next node
  // Find the correct next node chosen by the paths
  for(i = 0; i < num_next; i++)
    if(bases[i] == greatest_nuc)
      return_step(indices[i], GRPHWLK_USEPATH);

  // Should be impossible to reach here...
  assert(0);
}

#undef return_step



static void _graph_walker_pickup_counter_paths(GraphWalker *wlk,
                                               Nucleotide prev_nuc)
{
  const dBGraph *db_graph = wlk->db_graph;
  dBNode prev_nodes[4];
  Nucleotide prev_bases[4];
  size_t i, j, num_prev_nodes;
  Edges edges, prev_edge;
  Nucleotide next_base;
  Orientation backwards = !wlk->node.orient;

  // Remove edge to kmer we came from
  edges = db_node_get_edges(db_graph, 0, wlk->node.key);

  // Can slim down the number of nodes to look up if we can rule out
  // the node we just came from
  prev_nuc = dna_nuc_complement(prev_nuc);
  prev_edge = nuc_orient_to_edge(prev_nuc, backwards);

  // status("lost: %c:%i", dna_nuc_to_char(prev_nuc), backwards);

  // Some sanity checks
  assert(edges & prev_edge);
  assert(binary_kmers_are_equal(wlk->bkey, db_node_get_bkmer(db_graph, wlk->node.key)));

  num_prev_nodes = db_graph_next_nodes(db_graph, wlk->bkey,
                                       backwards, edges & ~prev_edge,
                                       prev_nodes, prev_bases);

  // If we have the ability, slim down nodes by those in this colour
  if(db_graph->node_in_cols != NULL) {
    for(i = j = 0; i < num_prev_nodes; i++) {
      if(db_node_has_col(db_graph, prev_nodes[i].key, wlk->ctxcol)) {
        prev_nodes[j] = prev_nodes[i];
        prev_bases[j] = prev_bases[i];
        j++;
      }
    }
    num_prev_nodes = j;
  }

  next_base = binary_kmer_last_nuc(wlk->bkmer);

  // Reverse orientation, pick up paths
  for(i = 0; i < num_prev_nodes; i++)
    pickup_paths(wlk, db_node_reverse(prev_nodes[i]), true, next_base);
}


static void _graph_traverse_force_jump(GraphWalker *wlk, hkey_t hkey,
                                       BinaryKmer bkmer, boolean is_fork)
{
  assert(hkey != HASH_NOT_FOUND);

  // #ifdef DEBUG_WALKER
  //   char str[MAX_KMER_SIZE+1];
  //   binary_kmer_to_str(bkmer, wlk->db_graph->kmer_size, str);
  //   printf("FORCE JUMP %s (fork:%s)\n", str, fork ? "yes" : "no");
  // #endif

  if(is_fork)
  {
    // We passed a fork - take all paths that agree with said nucleotide and
    // haven't ended, also update ages
    Nucleotide base = binary_kmer_last_nuc(bkmer), pnuc;
    FollowPath *path;
    size_t i, j, npaths = wlk->paths.len;

    path_buf_ensure_capacity(&wlk->paths, wlk->paths.len + wlk->new_paths.len);

    // Check curr paths
    path_buf_reset(&wlk->paths);

    for(i = 0, j = 0; i < npaths; i++)
    {
      path = &wlk->paths.data[i];
      pnuc = cache_fetch(path, path->pos);
      if(base == pnuc && path->pos+1 < path->len) {
        path->pos++;
        cache_update(path);
        wlk->paths.data[j++] = *path;
      }
    }

    // New paths -> curr paths
    for(i = 0; i < wlk->new_paths.len; i++)
    {
      path = &wlk->new_paths.data[i];
      pnuc = cache_fetch(path, path->pos);
      if(base == pnuc && path->pos+1 < path->len) {
        path->pos++;
        cache_update(path);
        wlk->paths.data[j++] = *path;
      }
    }

    wlk->paths.len = j;
    path_buf_reset(&wlk->new_paths);

    // Counter paths
    for(i = 0, j = 0; i < wlk->cntr_paths.len; i++)
    {
      path = &wlk->cntr_paths.data[i];
      pnuc = cache_fetch(path, path->pos);
      if(base == pnuc && path->pos+1 < path->len) {
        path->pos++;
        cache_update(path);
        wlk->cntr_paths.data[j++] = *path;
      }
    }

    wlk->cntr_paths.len = j;

    // Statistics
    wlk->fork_count++;
  }
  else if(wlk->new_paths.len > 0)
  {
    // Merge in (previously new) paths
    // New paths -> curr paths (no filtering required)
    path_buf_ensure_capacity(&wlk->paths, wlk->paths.len + wlk->new_paths.len);
    size_t mem = wlk->new_paths.len * sizeof(FollowPath);
    memcpy(wlk->paths.data+wlk->paths.len, wlk->new_paths.data, mem);
    wlk->paths.len += wlk->new_paths.len;
    path_buf_reset(&wlk->new_paths);
  }

  // Update GraphWalker position
  wlk->node.key = hkey;
  wlk->bkmer = bkmer;
  wlk->bkey = db_node_get_bkmer(wlk->db_graph, hkey);
  wlk->node.orient = db_node_get_orientation(wlk->bkmer, wlk->bkey);

  // Pick up new paths
  pickup_paths(wlk, wlk->node, false, 0);
}

// Jump to a new node within the current sample supernode
// (can actually be any node up until the end of the current supernode)
// If fork is true, node is the result of taking a fork -> slim down paths
// prev is array of nodes with edges to the node we are moving to
void graph_walker_jump_snode_end(GraphWalker *wlk, hkey_t hkey, BinaryKmer bkmer)
{
  // This is just a sanity test
  Edges edges = db_node_get_edges(wlk->db_graph, 0, hkey);
  BinaryKmer bkey = db_node_get_bkmer(wlk->db_graph, hkey);
  Orientation orient = db_node_get_orientation(bkmer, bkey);
  assert(edges_get_indegree(edges, orient) <= 1);

  // Now do the work
  _graph_traverse_force_jump(wlk, hkey, bkmer, false);
}

void graph_traverse_force(GraphWalker *wlk, hkey_t hkey, Nucleotide base,
                          boolean is_fork)
{
  assert(hkey != HASH_NOT_FOUND);
  BinaryKmer bkmer;
  const size_t kmer_size = wlk->db_graph->kmer_size;
  Nucleotide lost_nuc = binary_kmer_first_nuc(wlk->bkmer, kmer_size);
  bkmer = binary_kmer_left_shift_add(wlk->bkmer, kmer_size, base);

  // char tmp0[MAX_KMER_SIZE+1], tmp1[MAX_KMER_SIZE+1];
  // Orientation or0, or1;
  // binary_kmer_to_str(wlk->bkmer, kmer_size, tmp0);
  // binary_kmer_to_str(bkmer, kmer_size, tmp1);
  // status("%s -> %s lost: %c", tmp0, tmp1, dna_nuc_to_char(lost_nuc));
  // binary_kmer_to_str(wlk->bkey, kmer_size, tmp0);
  // binary_kmer_to_str(db_node_get_bkmer(wlk->db_graph,hkey), kmer_size, tmp1);
  // or0 = wlk->node.orient;
  // or1 = db_node_get_orientation(bkmer, db_node_get_bkmer(wlk->db_graph,hkey));
  // status("%s:%i -> %s:%i", tmp0, or0, tmp1, or1);

  _graph_traverse_force_jump(wlk, hkey, bkmer, is_fork);
  _graph_walker_pickup_counter_paths(wlk, lost_nuc);
}

boolean graph_traverse_nodes(GraphWalker *wlk, size_t num_next,
                             const dBNode nodes[4], const Nucleotide bases[4])
{
  wlk->last_step = graph_walker_choose(wlk, num_next, nodes, bases);
  int idx = wlk->last_step.idx;
  if(idx == -1) return false;
  graph_traverse_force(wlk, nodes[idx].key, bases[idx],
                       graphstep_is_fork(wlk->last_step));
  return true;
}

// return 1 on success, 0 otherwise
boolean graph_traverse(GraphWalker *wlk)
{
  const dBGraph *db_graph = wlk->db_graph;
  Edges edges = db_node_get_edges(db_graph, 0, wlk->node.key);

  dBNode nodes[4];
  Nucleotide bases[4];
  size_t num_next;

  num_next = db_graph_next_nodes(db_graph, wlk->bkey, wlk->node.orient, edges,
                                 nodes, bases);

  return graph_traverse_nodes(wlk, num_next, nodes, bases);
}


//
// Force traversal along an array of nodes (`priming' a GraphWalker)
//

// Fast traverse - avoid a bkmer_revcmp
static inline void graph_walker_fast(GraphWalker *wlk, const dBNode prev_node,
                                     const dBNode next_node, boolean is_fork)
{
  const size_t kmer_size = wlk->db_graph->kmer_size;
  BinaryKmer bkmer = db_node_get_bkmer(wlk->db_graph, next_node.key);
  Nucleotide nuc;

  // Only one path between two nodes
  if(db_nodes_match(wlk->node, prev_node)) {
    nuc = db_node_get_last_nuc(bkmer, next_node.orient, kmer_size);
    graph_traverse_force(wlk, next_node.key, nuc, is_fork);
  }
  else {
    // jumping to the end of a supernode
    bkmer = db_node_oriented_bkmer(bkmer, next_node.orient, kmer_size);
    graph_walker_jump_snode_end(wlk, next_node.key, bkmer);
  }

  // char tmpbkmer[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(wlk->bkmer, wlk->db_graph->kmer_size, tmpbkmer);
  // printf("  forced: %s\n", tmpbkmer);
}

// Fast traversal of a list of nodes using the supplied GraphWalker
// Only visits nodes deemed informative + last node
// Must have previously initialised or walked to the prior node,
// using: graph_walker_init, graph_traverse_force, graph_walker_jump_snode_end,
// graph_traverse or graph_traverse_nodes
// i.e. wlk->node is a node adjacent to arr[0]
void graph_walker_fast_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                boolean forward)
{
  if(n == 0) return;
  // Only one colour should be loaded
  assert(wlk->db_graph->num_of_cols == 1);

  size_t i;
  boolean infork[3] = {false, false, false}, outfork[3] = {false, false, false};
  Edges edges;
  dBNode nodes[3];

  edges = db_node_get_edges(wlk->db_graph, 0, wlk->node.key);
  outfork[0] = edges_get_outdegree(edges, wlk->node.orient) > 1;

  nodes[0] = wlk->node;
  nodes[1] = forward ? arr[0] : db_node_reverse(arr[n-1]);

  edges = db_node_get_edges(wlk->db_graph, 0, nodes[1].key);
  outfork[1] = edges_get_outdegree(edges, nodes[1].orient) > 1;
  infork[1] = edges_get_indegree(edges, nodes[1].orient) > 1;

  for(i = 0; i+1 < n; i++)
  {
    // printf("i: %zu %zu:%i\n", i, (size_t)nodes[1].key, (int)nodes[1].orient);
    nodes[2] = forward ? arr[i+1] : db_node_reverse(arr[n-i-2]);

    edges = db_node_get_edges(wlk->db_graph, 0, nodes[2].key);
    outfork[2] = edges_get_outdegree(edges, nodes[2].orient) > 1;
    infork[2] = edges_get_indegree(edges, nodes[2].orient) > 1;

    // Traverse nodes[i] if:
    // - previous node had out-degree > 1 (update/drop paths)
    // - current node has in-degree > 1 (pick up counter-paths + merge in new paths)
    // - next node has in-degree > 1 (pickup paths)
    if(outfork[0] || infork[1] || infork[2]) {
      graph_walker_fast(wlk, nodes[0], nodes[1], outfork[0]);
    }

    // Rotate edges, nodes
    infork[0] = infork[1]; infork[1] = infork[2];
    outfork[0] = outfork[1]; outfork[1] = outfork[2];
    nodes[0] = nodes[1]; nodes[1] = nodes[2];
  }

  // Traverse last node
  graph_walker_fast(wlk, nodes[0], nodes[1], outfork[0]);
}

// Traversal of every node in a list of nodes using the supplied GraphWalker
// Visits each node specifed
void graph_walker_slow_traverse(GraphWalker *wlk, const dBNode *arr, size_t n,
                                boolean forward)
{
  Edges edges;
  boolean is_fork;
  BinaryKmer bkmer;
  dBNode next;
  Nucleotide nuc;
  size_t i;
  const dBGraph *db_graph = wlk->db_graph;
  const size_t kmer_size = db_graph->kmer_size;

  for(i = 0; i < n; i++) {
    edges = db_node_get_edges(db_graph, 0, wlk->node.key);
    is_fork = edges_get_outdegree(edges, wlk->node.orient) > 1;
    next = forward ? arr[i] : db_node_reverse(arr[n-1-i]);
    bkmer = db_node_get_bkmer(db_graph, next.key);
    nuc = db_node_get_last_nuc(bkmer, next.orient, kmer_size);
    graph_traverse_force(wlk, next.key, nuc, is_fork);
  }
}
