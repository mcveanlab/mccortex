#include "global.h"
#include "graph_cache.h"
#include "binary_seq.h"
#include "supernode.h"

#include "sort_r/sort_r.h"

// Build and store paths through the graph
// Must build one path at a time
// Cannot update an old path
// (this could be overcome by making paths steps a linkedlist
//  i.e. adding uint32_t next_in_path field to step, but is not needed atm)
// Warning: not thread safe! Do not use the same GraphCache in more than one
//          thread at the same time.

void graph_cache_alloc(GraphCache *cache, const dBGraph *db_graph)
{
  db_node_buf_alloc(&cache->node_buf, 1024);
  cache_snode_buf_alloc(&cache->snode_buf, 1024);
  cache_step_buf_alloc(&cache->step_buf, 1024);
  cache_path_buf_alloc(&cache->path_buf, 1024);
  cache->snode_hash = kh_init(SnodeIdHash);
  cache->db_graph = db_graph;
}

void graph_cache_dealloc(GraphCache *cache)
{
  kh_destroy(SnodeIdHash, cache->snode_hash);
  db_node_buf_dealloc(&cache->node_buf);
  cache_snode_buf_dealloc(&cache->snode_buf);
  cache_step_buf_dealloc(&cache->step_buf);
  cache_path_buf_dealloc(&cache->path_buf);
}

// Encode supernode and orientation into a 32bit integer
#define _graph_cache_step_encode(step) (((step)->supernode << 1) | (step)->orient)

// Returns pathid
uint32_t graph_cache_new_path(GraphCache *cache)
{
  GCachePath path = {.first_step = cache->step_buf.len, .num_steps = 0};
  return cache_path_buf_add(&cache->path_buf, path);
}

// Create a supernode starting at node/or.  Store in snode.
// Ensure snode->nodes and snode->orients point to valid memory before passing
// Returns 0 on failure, otherwise snode->num_of_nodes
static inline void _create_supernode(GraphCache *cache, dBNode node,
                                     GCacheSnode *snode)
{
  const dBGraph *db_graph = cache->db_graph;
  ctx_assert(db_graph->num_edge_cols == 1);

  size_t first_node_id = cache->node_buf.len;
  db_node_buf_add(&cache->node_buf, node);
  supernode_extend(&cache->node_buf, 0, db_graph);
  size_t num_nodes = cache->node_buf.len - first_node_id;

  dBNode *nodes = graph_cache_node(cache, first_node_id);
  supernode_normalise(nodes, num_nodes, db_graph);

  // printf("Loaded supernode:\n  ");
  // db_nodes_print(nodes, num_nodes, db_graph, stdout);
  // printf("\n");

  BinaryKmer bkmer0, bkmer1;
  Edges union_edges;
  dBNode prev_nodes[4], next_nodes[4];
  Nucleotide prev_bases[4] = {0}, next_bases[4] = {0};
  uint8_t num_prev, num_next;

  dBNode first = db_node_reverse(nodes[0]);
  dBNode last = nodes[num_nodes-1];

  // prev nodes
  union_edges = db_node_get_edges_union(db_graph, first.key);
  bkmer0 = db_node_get_bkmer(db_graph, first.key);
  num_prev = db_graph_next_nodes(db_graph, bkmer0, first.orient, union_edges,
                                 prev_nodes, prev_bases);

  // next nodes
  union_edges = db_node_get_edges_union(db_graph, last.key);
  bkmer1 = db_node_get_bkmer(db_graph, last.key);
  num_next = db_graph_next_nodes(db_graph, bkmer1, last.orient, union_edges,
                                 next_nodes, next_bases);

  uint8_t prev_packed = binary_seq_pack_byte(prev_bases);
  uint8_t next_packed = binary_seq_pack_byte(next_bases);

  GCacheSnode tmp = {.first_node_id = first_node_id,
                        .num_nodes = num_nodes,
                        .first_step = UINT32_MAX,
                        .num_prev = num_prev,
                        .prev_bases = prev_packed,
                        .prev_nodes[0] = prev_nodes[0],
                        .prev_nodes[1] = prev_nodes[1],
                        .prev_nodes[2] = prev_nodes[2],
                        .prev_nodes[3] = prev_nodes[3],
                        .num_next = num_next,
                        .next_bases = next_packed,
                        .next_nodes[0] = next_nodes[0],
                        .next_nodes[1] = next_nodes[1],
                        .next_nodes[2] = next_nodes[2],
                        .next_nodes[3] = next_nodes[3]};

  memcpy(snode, &tmp, sizeof(GCacheSnode));
}

Orientation graph_cache_get_supernode_orient(const GraphCache *cache,
                                             const GCacheSnode *snode,
                                             dBNode first_node)
{
  dBNode *snode0 = graph_cache_node(cache, snode->first_node_id);
  return db_nodes_are_equal(*snode0, first_node) ? FORWARD : REVERSE;
}

// Get node from opposite end of the supernode
static inline dBNode _node_at_snode_end(GraphCache *cache,
                                        GCacheSnode *snode,
                                        dBNode node)
{
  dBNode *first_node = graph_cache_node(cache, snode->first_node_id);
  if(db_nodes_are_equal(*first_node, node))
    return db_node_reverse(first_node[snode->num_nodes-1]);
  else
    return *first_node;
}

// Returns stepid
uint32_t graph_cache_new_step(GraphCache *cache, dBNode node)
{
  // Get current path
  uint32_t pathid = cache->path_buf.len-1;
  GCachePath *path = graph_cache_path(cache, pathid);

  // Find or add supernode beginning with given node
  uint32_t snodeid;
  int hashret;
  khiter_t k = kh_put(SnodeIdHash, cache->snode_hash, node, &hashret);
  bool supernode_already_exists = (hashret == 0);

  if(supernode_already_exists) {
    snodeid = kh_value(cache->snode_hash, k);
  }
  else {
    // Create supernode
    GCacheSnode tmp_snode;
    _create_supernode(cache, node, &tmp_snode);
    snodeid = cache_snode_buf_add(&cache->snode_buf, tmp_snode);
    kh_value(cache->snode_hash, k) = snodeid;

    // Get node at other end
    dBNode end_node = _node_at_snode_end(cache, &tmp_snode, node);
    k = kh_put(SnodeIdHash, cache->snode_hash, end_node, &hashret);
    kh_value(cache->snode_hash, k) = snodeid;
  }

  GCacheSnode *snode = graph_cache_snode(cache, snodeid);

  // Get orient
  Orientation snode_orient = graph_cache_get_supernode_orient(cache, snode, node);

  // New step
  GCacheStep next = {.orient = snode_orient, .supernode = snodeid,
                    .pathid = pathid, .next_step = snode->first_step};
  uint32_t stepid = cache_step_buf_add(&cache->step_buf, next);

  // Add link from prev supernode step
  snode->first_step =  stepid;

  path->num_steps++;
  return stepid;
}

void graph_cache_reset(GraphCache *cache)
{
  kh_clear(SnodeIdHash, cache->snode_hash);
  db_node_buf_reset(&cache->node_buf);
  cache_snode_buf_reset(&cache->snode_buf);
  cache_step_buf_reset(&cache->step_buf);
  cache_path_buf_reset(&cache->path_buf);
}

int graph_cache_steps_cmp(const GCacheStep *a, const GCacheStep *b,
                          const GraphCache *cache)
{
  // Sort by supernodes up to (and including) this point
  const GCachePath *apath = graph_cache_path(cache, a->pathid);
  const GCachePath *bpath = graph_cache_path(cache, b->pathid);

  const GCacheStep *step0 = graph_cache_step(cache, apath->first_step);
  const GCacheStep *step1 = graph_cache_step(cache, bpath->first_step);

  uint32_t len0 = a - step0 + 1, len1 = b - step1 + 1, minlen = MIN2(len0, len1);

  // compare path steps
  uint32_t word0, word1;
  size_t i;

  for(i = 0; i < minlen; i++) {
    word0 = _graph_cache_step_encode(step0);
    word1 = _graph_cache_step_encode(step1);
    if(word0 < word1) return -1;
    if(word0 > word1) return 1;
  }

  return (long)len0 - len1;
}

static inline int stepptr_cmp(const void *aa, const void *bb, void *arg)
{
  const GCacheStep *const*a = (const GCacheStep *const*)aa;
  const GCacheStep *const*b = (const GCacheStep *const*)bb;
  return graph_cache_steps_cmp(*a, *b, (const GraphCache *)arg);
}

void graph_cache_stepptrs_qsort(GraphCache *cache, GCacheStep **list, size_t n)
{
  sort_r(list, n, sizeof(GCacheStep*), stepptr_cmp, cache);
}

static inline int pathids_cmp(const void *aa, const void *bb, void *arg)
{
  uint32_t a = *(const uint32_t *)aa, b = *(const uint32_t *)bb;
  const GraphCache *cache = (const GraphCache*)arg;
  const GCachePath *patha = graph_cache_path(cache, a);
  const GCachePath *pathb = graph_cache_path(cache, b);
  const GCacheStep *stepa = graph_cache_path_last_step(cache, patha);
  const GCacheStep *stepb = graph_cache_path_last_step(cache, pathb);
  return graph_cache_steps_cmp(stepa, stepb, cache);
}

bool graph_cache_pathids_are_equal(GraphCache *cache,
                                   uint32_t pathid0, uint32_t pathid1)
{
  return pathids_cmp(&pathid0, &pathid1, cache) == 0;
}

// Get all nodes in a single step (supernode with orientation)
// Adds to the end of the node buffer (does not reset it)
void graph_cache_snode_fetch_nodes(const GraphCache *cache,
                                   const GCacheSnode *snode,
                                   Orientation orient,
                                   dBNodeBuffer *nbuf)
{
  const dBNode *nodes = graph_cache_node(cache, snode->first_node_id);
  const dBNode *end = nodes + snode->num_nodes;
  dBNode *into;

  db_node_buf_ensure_capacity(nbuf, nbuf->len + snode->num_nodes);

  if(orient == FORWARD) {
    memcpy(nbuf->data + nbuf->len, nodes, sizeof(dBNode) * snode->num_nodes);
  }
  else {
    for(into = nbuf->data+nbuf->len+snode->num_nodes-1; nodes < end; into--, nodes++)
      *into = db_node_reverse(*nodes);
  }

  nbuf->len += snode->num_nodes;
}

// Get all nodes in a path up to, but not including the given step
// Adds to the end of the node buffer (does not reset it)
void graph_cache_step_fetch_nodes(const GraphCache *cache,
                                  const GCacheStep *endstep,
                                  dBNodeBuffer *nbuf)
{
  const GCachePath *path = graph_cache_path(cache, endstep->pathid);
  const GCacheStep *step = graph_cache_step(cache, path->first_step);
  const GCacheSnode *snode;

  // Loop over steps, load nodes
  for(; step < endstep; step++) {
    snode = graph_cache_snode(cache, step->supernode);
    graph_cache_snode_fetch_nodes(cache, snode, step->orient, nbuf);
  }
}

void graph_cache_path_fetch_nodes(const GraphCache *cache,
                                  const GCachePath *path,
                                  dBNodeBuffer *nbuf)
{
  const GCacheStep *step = graph_cache_path_last_step(cache, path);
  graph_cache_step_fetch_nodes(cache, step+1, nbuf);
}

// Get previous step or NULL if first step
static const GCacheStep* _fetch_prev_step(GraphCache *cache, GCacheStep *step)
{
  GCachePath *path = graph_cache_path(cache, step->pathid);
  GCacheStep *step0 = graph_cache_step(cache, path->first_step);
  return (step0 == step ? NULL : --step);
}

// Get first step in path
static const GCacheStep* _fetch_first_step(GraphCache *cache, GCacheStep *step)
{
  GCachePath *path = graph_cache_path(cache, step->pathid);
  return graph_cache_step(cache, path->first_step);
}

// Looks like 3p flank if steps don't have the same n-1 supernode
bool graph_cache_is_3p_flank(GraphCache *cache,
                             GCacheStep **steps, size_t num_steps)
{
  if(num_steps <= 1) return false;

  const GCacheStep *step0, *step1;
  size_t i;
  uint32_t word0, word1;

  // 1. Check first step differs
  const GCacheStep *first1, *first0 = _fetch_first_step(cache, steps[0]);
  word0 = _graph_cache_step_encode(first0);

  for(i = 1; i < num_steps; i++) {
    first1 = _fetch_first_step(cache, steps[1]);
    word1 = _graph_cache_step_encode(first1);
    if(word0 != word1) break;
  }

  //   Not a bubble if all paths start with the same node
  if(i == num_steps) return false;

  // 2. Check second last step differs (supernode before 3p flank)
  step0 = _fetch_prev_step(cache, steps[0]);

  // If first step is null, we only have to find one that is not-null
  if(step0 == NULL) {
    for(i = 1; i < num_steps; i++)
      if(_fetch_prev_step(cache, steps[i]) != NULL)
        return true;
  }
  else {
    // Have to find a prev step that is null or has diff supernode
    word0 = _graph_cache_step_encode(step0);
    for(i = 1; i < num_steps; i++) {
      step1 = _fetch_prev_step(cache, steps[i]);
      word1 = _graph_cache_step_encode(step1);
      if(step1 == NULL || word0 != word1)
        return true;
    }
  }

  return false;
}

// Remove duplicate paths
size_t graph_cache_remove_dupes(GraphCache *cache,
                                GCacheStep **steps, size_t num_steps)
{
  size_t i, j;

  if(num_steps <= 1) return num_steps;

  graph_cache_stepptrs_qsort(cache, steps, num_steps);

  for(i = j = 0; i+1 < num_steps; i++) {
    if(graph_cache_steps_cmp(steps[i], steps[i+1], cache) != 0)
      steps[j++] = steps[i];
  }

  steps[j++] = steps[i];
  return j;
}

// Returns true if all nodes in supernode have given colour
bool graph_cache_snode_has_colour(const GraphCache *cache,
                                  const GCacheSnode *snode,
                                  size_t colour)
{
  const dBGraph *db_graph = cache->db_graph;
  const dBNode *node = graph_cache_node(cache, snode->first_node_id), *end;
  for(end = node + snode->num_nodes; node < end; node++) {
    if(!db_node_has_col(db_graph, node->key, colour))
      return false;
  }
  return true;
}

// Returns true if all nodes in path have given colour
bool graph_cache_step_has_colour(const GraphCache *cache,
                                 const GCacheStep *endstep,
                                 size_t colour)
{
  const GCacheStep *step = graph_cache_step(cache, endstep->pathid);
  const GCacheSnode *snode;

  for(; step <= endstep; step++) {
    snode = graph_cache_snode(cache, step->supernode);
    if(!graph_cache_snode_has_colour(cache, snode, colour)) return false;
  }

  return true;
}

// Returns NULL if not found
GCacheSnode* graph_cache_find_snode(GraphCache *cache, dBNode node)
{
  uint32_t snodeid;
  khiter_t k = kh_get(SnodeIdHash, cache->snode_hash, node);
  bool is_missing = (k == kh_end(cache->snode_hash));

  if(!is_missing) {
    snodeid = kh_value(cache->snode_hash, k);
    return graph_cache_snode(cache, snodeid);
  }

  return NULL;
}
