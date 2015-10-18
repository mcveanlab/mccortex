#ifndef GRAPH_CACHE_H_
#define GRAPH_CACHE_H_

#include "htslib/khash.h"
#include "db_node.h"

// Build and store paths through the graph
// Must build one path at a time
// Cannot update an old path
// (this could be overcome by making paths steps a linkedlist
//  i.e. adding uint32_t next_in_path field to step, but is not needed atm)
// Warning: not thread safe! Do not use the same GraphCache in more than one
//          thread at the same time.

// typedef uint32_t GCUnitigId, GCStepId, GCPathId;

typedef struct
{
  const size_t first_node_id;
  uint32_t num_nodes;
  uint32_t stepid; // linked list of steps through this unitig

  // This is to speed up traversal of subsequent colours
  const dBNode prev_nodes[4], next_nodes[4];
  const uint8_t prev_bases, next_bases; // bases packed into 2bits per base
  const uint8_t num_prev:4, num_next:4;
} GCacheUnitig;

// A step is a unitig with an orientation
typedef struct
{
  const uint32_t unitigid:31, orient:1;
  const uint32_t pathid; // path that this step belongs to
  uint32_t next_step; // linked list of steps that visit this unitig
} GCacheStep;

typedef struct
{
  const uint32_t first_step;
  uint32_t num_steps;
} GCachePath;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(cache_unitig_buf, GCacheUnitigBuffer, GCacheUnitig);
madcrow_buffer(cache_step_buf,   GCacheStepBuffer,   GCacheStep);
madcrow_buffer(cache_path_buf,   GCachePathBuffer,   GCachePath);

KHASH_INIT(Node2Unitig, dBNode, uint32_t, 1, db_node_hash, db_nodes_are_equal)

typedef struct
{
  dBNodeBuffer       node_buf;
  GCacheUnitigBuffer unitig_buf;
  GCacheStepBuffer   step_buf;
  GCachePathBuffer   path_buf;

  // hash map dBNode->uint32_t (unitig_id)
  khash_t(Node2Unitig) *node2unitig;

  const dBGraph *db_graph;
} GraphCache;

void graph_cache_alloc(GraphCache *cache, const dBGraph *db_graph);
void graph_cache_dealloc(GraphCache *cache);
void graph_cache_reset(GraphCache *cache);

#define graph_cache_node(cache,nodeid) (&(cache)->node_buf.b[nodeid])
#define graph_cache_unitig(cache,unitigid) (&(cache)->unitig_buf.b[unitigid])
#define graph_cache_step(cache,stepid) (&(cache)->step_buf.b[stepid])
#define graph_cache_path(cache,pathid) (&(cache)->path_buf.b[pathid])

#define graph_cache_node_id(cache,node) ((node) - (cache)->node_buf.b)
#define graph_cache_unitig_id(cache,unitig) ((unitig) - (cache)->unitig_buf.b)
#define graph_cache_step_id(cache,step) ((step) - (cache)->step_buf.b)
#define graph_cache_path_id(cache,path) ((path) - (cache)->path_buf.b)

#define graph_cache_num_nodes(cache) ((cache)->node_buf.len)
#define graph_cache_num_unitigs(cache) ((cache)->unitig_buf.len)
#define graph_cache_num_steps(cache) ((cache)->step_buf.len)
#define graph_cache_num_paths(cache) ((cache)->path_buf.len)

//
// BEWARE: we resize buffers if needed each time we add a new unitig/step/path
// Calling graph_cache_new_path() invalidates all pointers to paths in the cache
// except for the pointer returned.
// The same is true from graph_cache_new_step() and GCacheStep pointers.
//

// Returns new path
const GCachePath* graph_cache_new_path(GraphCache *cache);

// Returns new step
const GCacheStep* graph_cache_new_step(GraphCache *cache, dBNode node);

// Returns NULL if not found
GCacheUnitig* graph_cache_find_unitig(GraphCache *cache, dBNode node);

//
// Paths
//

static inline const GCacheStep* gc_path_first_step(const GraphCache *cache,
                                                   const GCachePath *path)
{
  return graph_cache_step(cache, path->first_step);
}

static inline const GCacheStep* gc_path_last_step(const GraphCache *cache,
                                                  const GCachePath *path)
{
  return graph_cache_step(cache, path->first_step + path->num_steps - 1);
}

// Get number of kmers in a path
size_t gc_path_get_nkmers(const GraphCache *cache, const GCachePath *path);

// Get all nodes in a path
// Adds to the end of the node buffer (does not reset it)
void gc_path_fetch_nodes(const GraphCache *cache,
                         const GCachePath *path, size_t num_steps,
                         dBNodeBuffer *nbuf);

//
// Unitigs
//

static inline const dBNode* gc_unitig_get_nodes(const GraphCache *cache,
                                                const GCacheUnitig *u)
{
  return graph_cache_node(cache, u->first_node_id);
}

static inline dBNode gc_unitig_first_node(const GraphCache *cache,
                                          const GCacheUnitig *u)
{
  return *graph_cache_node(cache, u->first_node_id);
}

static inline dBNode gc_unitig_last_node(const GraphCache *cache,
                                         const GCacheUnitig *u)
{
  return *graph_cache_node(cache, u->first_node_id + u->num_nodes - 1);
}

Orientation gc_unitig_get_orient(const GraphCache *cache,
                                 const GCacheUnitig *unitig,
                                 dBNode first_node);

// Get all nodes in a single step (unitig with orientation)
// Adds to the end of the node buffer (does not reset it)
void gc_unitig_fetch_nodes(const GraphCache *cache,
                           const GCacheUnitig *unitig,
                           Orientation orient,
                           dBNodeBuffer *nbuf);

//
// Steps
//

static inline const GCacheUnitig* gc_step_get_unitig(const GraphCache *cache,
                                                     const GCacheStep *s)
{
  return graph_cache_unitig(cache, s->unitigid);
}

static inline const GCachePath* gc_step_get_path(const GraphCache *cache,
                                                 const GCacheStep *s)
{
  return graph_cache_path(cache, s->pathid);
}

static inline const GCacheStep* gc_step_get_next(const GraphCache *cache,
                                                 const GCacheStep *s)
{
  return s->next_step == UINT32_MAX ? NULL : graph_cache_step(cache, s->next_step);
}

// Get previous step or NULL if first step
static inline const GCacheStep* gc_step_get_prev(const GraphCache *cache,
                                                 const GCacheStep *step)
{
  const GCacheStep *step0 = gc_path_first_step(cache, gc_step_get_path(cache, step));
  return (step0 == step ? NULL : step-1);
}

// Encode unitig and orientation into a 32bit integer
static inline uint32_t gc_step_encode_uint32(const GCacheStep *step)
{
  return (((uint32_t)step->unitigid << 1) | step->orient);
}

static inline size_t gc_step_get_nkmers(const GraphCache *g, const GCacheStep *s) {
  return graph_cache_unitig(g, s->unitigid)->num_nodes;
}

// Get all nodes in a path up to, but not including the given step
// Adds to the end of the node buffer (does not reset it)
void gc_step_fetch_nodes(const GraphCache *cache,
                         const GCacheStep *end_step,
                         dBNodeBuffer *nbuf);

//
// Sorting
//

int graph_cache_pathids_cmp(const void *aa, const void *bb, void *arg);

int graph_cache_steps_cmp(const GCacheStep *a, const GCacheStep *b,
                          const GraphCache *cache);

void graph_cache_steps_qsort(GraphCache *cache, GCacheStep **list, size_t n);

bool graph_cache_pathids_are_equal(GraphCache *cache,
                                   uint32_t pathid0, uint32_t pathid1);

//
// Variant calling
//

// Looks like 3p flank if steps don't have the same n-1 unitig
bool graph_cache_is_3p_flank(GraphCache *cache,
                             GCacheStep ** steps, size_t num_steps);

// Remove duplicate paths
size_t graph_cache_remove_dupes(GraphCache *cache,
                                GCacheStep **steps, size_t num_steps);

// Returns true if all nodes in unitig have given colour
bool graph_cache_unitig_has_colour(const GraphCache *cache,
                                  const GCacheUnitig *unitig,
                                  size_t colour);

// Returns true if all nodes in path have given colour
bool graph_cache_step_has_colour(const GraphCache *cache,
                                 const GCacheStep *endstep,
                                 size_t colour);

#endif /* GRAPH_CACHE_H_ */
