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
  cache_unitig_buf_alloc(&cache->unitig_buf, 1024);
  cache_step_buf_alloc(&cache->step_buf, 1024);
  cache_path_buf_alloc(&cache->path_buf, 1024);
  cache->node2unitig = kh_init(Node2Unitig);
  cache->db_graph = db_graph;
}

void graph_cache_dealloc(GraphCache *cache)
{
  kh_destroy(Node2Unitig, cache->node2unitig);
  db_node_buf_dealloc(&cache->node_buf);
  cache_unitig_buf_dealloc(&cache->unitig_buf);
  cache_step_buf_dealloc(&cache->step_buf);
  cache_path_buf_dealloc(&cache->path_buf);
}

void graph_cache_reset(GraphCache *cache)
{
  kh_clear(Node2Unitig, cache->node2unitig);
  db_node_buf_reset(&cache->node_buf);
  cache_unitig_buf_reset(&cache->unitig_buf);
  cache_step_buf_reset(&cache->step_buf);
  cache_path_buf_reset(&cache->path_buf);
}

// Returns pathid
const GCachePath* graph_cache_new_path(GraphCache *cache)
{
  GCachePath path = {.first_step = cache->step_buf.len, .num_steps = 0};
  uint32_t pid = cache_path_buf_add(&cache->path_buf, path);
  return cache->path_buf.b + pid;
}

// Create a unitig starting at node/or.  Store in unitig.
// Ensure unitig->nodes and unitig->orients point to valid memory before passing
// Returns 0 on failure, otherwise unitig->num_of_nodes
static inline void gc_create_unitig(GraphCache *cache, dBNode node,
                                    GCacheUnitig *unitig)
{
  const dBGraph *db_graph = cache->db_graph;
  ctx_assert(db_graph->num_edge_cols == 1);

  size_t first_node_id = cache->node_buf.len;
  db_node_buf_add(&cache->node_buf, node);
  supernode_extend(&cache->node_buf, 0, db_graph);
  size_t num_nodes = cache->node_buf.len - first_node_id;

  dBNode *nodes = graph_cache_node(cache, first_node_id);
  supernode_normalise(nodes, num_nodes, db_graph);

  // printf("Loaded unitig:\n  ");
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

  GCacheUnitig tmp = {.first_node_id = first_node_id,
                     .num_nodes = num_nodes,
                     .stepid = UINT32_MAX,
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

  memcpy(unitig, &tmp, sizeof(GCacheUnitig));
}

// Get node from opposite end of the unitig
static inline dBNode get_node_at_unitig_end(GraphCache *cache,
                                            GCacheUnitig *unitig,
                                            dBNode node)
{
  dBNode *first_node = graph_cache_node(cache, unitig->first_node_id);
  if(db_nodes_are_equal(*first_node, node))
    return db_node_reverse(first_node[unitig->num_nodes-1]);
  else
    return *first_node;
}

// Returns new step
const GCacheStep* graph_cache_new_step(GraphCache *cache, dBNode node)
{
  // Get current path
  uint32_t pathid = cache->path_buf.len-1;
  GCachePath *path = graph_cache_path(cache, pathid);

  // Find or add unitig beginning with given node
  uint32_t unitigid;
  int hashret;
  khiter_t k = kh_put(Node2Unitig, cache->node2unitig, node, &hashret);
  bool unitig_already_exists = (hashret == 0);

  if(unitig_already_exists) {
    unitigid = kh_value(cache->node2unitig, k);
  }
  else {
    // Create unitig
    GCacheUnitig tmp_unitig;
    gc_create_unitig(cache, node, &tmp_unitig);
    unitigid = cache_unitig_buf_add(&cache->unitig_buf, tmp_unitig);
    kh_value(cache->node2unitig, k) = unitigid;

    // Get node at other end
    dBNode end_node = get_node_at_unitig_end(cache, &tmp_unitig, node);
    k = kh_put(Node2Unitig, cache->node2unitig, end_node, &hashret);
    kh_value(cache->node2unitig, k) = unitigid;
  }

  GCacheUnitig *unitig = graph_cache_unitig(cache, unitigid);

  // Get orient
  Orientation unitig_orient = gc_unitig_get_orient(cache, unitig, node);

  // New step
  GCacheStep next = {.orient = unitig_orient, .unitigid = unitigid,
                     .pathid = pathid, .next_step = unitig->stepid};
  uint32_t stepid = cache_step_buf_push(&cache->step_buf, &next, 1);

  // Add link from prev unitig step
  unitig->stepid = stepid;

  path->num_steps++;
  return cache->step_buf.b + stepid;
}

// Returns NULL if not found
GCacheUnitig* graph_cache_find_unitig(GraphCache *cache, dBNode node)
{
  uint32_t unitigid;
  khiter_t k = kh_get(Node2Unitig, cache->node2unitig, node);
  bool is_missing = (k == kh_end(cache->node2unitig));

  if(!is_missing) {
    unitigid = kh_value(cache->node2unitig, k);
    return graph_cache_unitig(cache, unitigid);
  }

  return NULL;
}


//
// Paths
//

size_t gc_path_get_nkmers(const GraphCache *cache, const GCachePath *path)
{
  size_t n = 0;
  const GCacheStep *endstep, *step = graph_cache_step(cache, path->first_step);
  for(endstep = step + path->num_steps; step < endstep; step++)
    n += gc_step_get_nkmers(cache, step);
  return n;
}

void gc_path_fetch_nodes(const GraphCache *cache,
                         const GCachePath *path, size_t num_steps,
                         dBNodeBuffer *nbuf)
{
  ctx_assert(num_steps <= path->num_steps);

  const GCacheStep *endstep, *step = graph_cache_step(cache, path->first_step);
  const GCacheUnitig *unitig;

  for(endstep = step + num_steps; step < endstep; step++) {
    unitig = graph_cache_unitig(cache, step->unitigid);
    gc_unitig_fetch_nodes(cache, unitig, step->orient, nbuf);
  }
}

//
// Unitigs
//

Orientation gc_unitig_get_orient(const GraphCache *cache,
                                 const GCacheUnitig *unitig,
                                 dBNode first_node)
{
  dBNode *unitig0 = graph_cache_node(cache, unitig->first_node_id);
  return db_nodes_are_equal(*unitig0, first_node) ? FORWARD : REVERSE;
}

// Get all nodes in a single step (unitig with orientation)
// Adds to the end of the node buffer (does not reset it)
void gc_unitig_fetch_nodes(const GraphCache *cache,
                           const GCacheUnitig *unitig,
                           Orientation orient,
                           dBNodeBuffer *nbuf)
{
  const dBNode *nodes = graph_cache_node(cache, unitig->first_node_id);
  const dBNode *end = nodes + unitig->num_nodes;
  dBNode *into;

  db_node_buf_capacity(nbuf, nbuf->len + unitig->num_nodes);

  if(orient == FORWARD) {
    memcpy(nbuf->b + nbuf->len, nodes, sizeof(dBNode) * unitig->num_nodes);
  }
  else {
    for(into = nbuf->b+nbuf->len+unitig->num_nodes-1; nodes < end; into--, nodes++)
      *into = db_node_reverse(*nodes);
  }

  nbuf->len += unitig->num_nodes;
}

//
// Steps
//

// Get all nodes in a path up to, but not including the given step
// Adds to the end of the node buffer (does not reset it)
void gc_step_fetch_nodes(const GraphCache *cache,
                         const GCacheStep *endstep,
                         dBNodeBuffer *nbuf)
{
  const GCachePath *path = graph_cache_path(cache, endstep->pathid);
  const GCacheStep *step0 = gc_path_first_step(cache, path);
  gc_path_fetch_nodes(cache, path, endstep-step0, nbuf);
}

//
// Sorting
//

// Sort by unitigs up to (and including) this point
int graph_cache_steps_cmp(const GCacheStep *a, const GCacheStep *b,
                          const GraphCache *cache)
{
  const GCachePath *path0, *path1;
  const GCacheStep *step0, *step1, *endstep0;

  path0 = graph_cache_path(cache, a->pathid);
  path1 = graph_cache_path(cache, b->pathid);

  step0 = graph_cache_step(cache, path0->first_step);
  step1 = graph_cache_step(cache, path1->first_step);

  // compare path steps
  uint32_t len0 = a - step0 + 1, len1 = b - step1 + 1, minlen = MIN2(len0, len1);
  uint32_t word0, word1;

  for(endstep0 = step0+minlen; step0 < endstep0; step0++, step1++) {
    word0 = gc_step_encode_uint32(step0);
    word1 = gc_step_encode_uint32(step1);
    if(word0 < word1) return -1;
    if(word0 > word1) return 1;
  }

  return cmp(len0, len1);
}

static inline int _steps_cmp(const void *aa, const void *bb, void *arg)
{
  const GCacheStep *const*a = (const GCacheStep *const*)aa;
  const GCacheStep *const*b = (const GCacheStep *const*)bb;
  return graph_cache_steps_cmp(*a, *b, (const GraphCache *)arg);
}

void graph_cache_steps_qsort(GraphCache *cache, GCacheStep **list, size_t n)
{
  sort_r(list, n, sizeof(GCacheStep*), _steps_cmp, cache);
}

int graph_cache_pathids_cmp(const void *aa, const void *bb, void *arg)
{
  uint32_t a = *(const uint32_t *)aa, b = *(const uint32_t *)bb;
  const GraphCache *cache = (const GraphCache*)arg;
  const GCachePath *patha = graph_cache_path(cache, a);
  const GCachePath *pathb = graph_cache_path(cache, b);
  if(patha->num_steps == 0 || pathb->num_steps == 0)
    return patha->num_steps - pathb->num_steps;
  const GCacheStep *stepa = gc_path_last_step(cache, patha);
  const GCacheStep *stepb = gc_path_last_step(cache, pathb);
  return graph_cache_steps_cmp(stepa, stepb, cache);
}

bool graph_cache_pathids_are_equal(GraphCache *cache,
                                   uint32_t pathid0, uint32_t pathid1)
{
  return graph_cache_pathids_cmp(&pathid0, &pathid1, cache) == 0;
}


//
// Variant calling
//

// Looks like 3p flank if steps don't have the same n-1 unitig
bool graph_cache_is_3p_flank(GraphCache *cache,
                             GCacheStep **steps, size_t num_steps)
{
  if(num_steps <= 1) return false;

  const GCacheStep *step0, *step1;
  size_t i;
  uint32_t word0, word1;

  // 1. Check first step differs
  const GCacheStep *first1, *first0;
  first0 = gc_path_first_step(cache, gc_step_get_path(cache, steps[0]));
  word0 = gc_step_encode_uint32(first0);

  for(i = 1; i < num_steps; i++) {
    first1 = gc_path_first_step(cache, gc_step_get_path(cache, steps[i]));
    word1 = gc_step_encode_uint32(first1);
    if(word0 != word1) break;
  }

  //   Not a bubble if all paths start with the same node
  if(i == num_steps) return false;

  // 2. Check second last step differs (unitig before 3p flank)
  step0 = gc_step_get_prev(cache, steps[0]);

  // If first step is null, we only have to find one that is not-null
  if(step0 == NULL) {
    for(i = 1; i < num_steps; i++)
      if(gc_step_get_prev(cache, steps[i]) != NULL)
        return true;
  }
  else {
    // Have to find a prev step that is null or has diff unitig
    word0 = gc_step_encode_uint32(step0);
    for(i = 1; i < num_steps; i++) {
      step1 = gc_step_get_prev(cache, steps[i]);
      if(step1 == NULL) return true;
      word1 = gc_step_encode_uint32(step1);
      if(word0 != word1) return true;
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

  graph_cache_steps_qsort(cache, steps, num_steps);

  for(i = j = 0; i+1 < num_steps; i++) {
    if(graph_cache_steps_cmp(steps[i], steps[i+1], cache) != 0)
      steps[j++] = steps[i];
  }

  steps[j++] = steps[i];
  return j;
}

// Returns true if all nodes in unitig have given colour
bool graph_cache_unitig_has_colour(const GraphCache *cache,
                                  const GCacheUnitig *unitig,
                                  size_t colour)
{
  const dBGraph *db_graph = cache->db_graph;
  const dBNode *node = graph_cache_node(cache, unitig->first_node_id), *end;
  for(end = node + unitig->num_nodes; node < end; node++) {
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
  const GCachePath *path = gc_step_get_path(cache, endstep);
  const GCacheStep *step = gc_path_first_step(cache, path);
  const GCacheUnitig *unitig;

  for(; step <= endstep; step++) {
    unitig = graph_cache_unitig(cache, step->unitigid);
    if(!graph_cache_unitig_has_colour(cache, unitig, colour)) return false;
  }

  return true;
}
