#ifndef GRAPH_CACHE_H_
#define GRAPH_CACHE_H_

#include "khash.h"
#include "db_node.h"

// Build and store paths through the graph
// Must build one path at a time
// Cannot update an old path
// (this could be overcome by making paths steps a linkedlist
//  i.e. adding uint32_t next_in_path field to step, but is not needed atm)
// Warning: not thread safe! Do not use the same GraphCache in more than one
//          thread at the same time.

typedef struct
{
  const size_t first_node_id;
  uint32_t num_nodes;
  uint32_t first_step; // linked list of steps through this supernode

  // This is to speed up traversal of subsequent colours
  const dBNode prev_nodes[4], next_nodes[4];
  const uint8_t prev_bases, next_bases; // bases packed into 2bits per base
  const uint8_t num_prev:4, num_next:4;
} GCacheSnode;

typedef struct
{
  const uint32_t orient:1, supernode:31;
  const uint32_t pathid; // path that this step belongs to
  uint32_t next_step; // linked list of steps through a single supernode
} GCacheStep;

typedef struct
{
  const uint32_t first_step;
  uint32_t num_steps;
} GCachePath;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(cache_snode_buf, GCacheSnodeBuffer, GCacheSnode);
madcrow_buffer(cache_step_buf,  GCacheStepBuffer,  GCacheStep);
madcrow_buffer(cache_path_buf,  GCachePathBuffer,  GCachePath);

#define db_node_hash(x) kh_int64_hash_func((x.key << 1) | x.orient)
KHASH_INIT(SnodeIdHash, dBNode, uint32_t, 1, db_node_hash, db_nodes_are_equal)

typedef struct
{
  dBNodeBuffer       node_buf;
  GCacheSnodeBuffer  snode_buf;
  GCacheStepBuffer   step_buf;
  GCachePathBuffer   path_buf;

  // hash hkey_t->uint32_t (supernode_id)
  khash_t(SnodeIdHash) *snode_hash;

  const dBGraph *db_graph;
} GraphCache;

void graph_cache_alloc(GraphCache *cache, const dBGraph *db_graph);
void graph_cache_dealloc(GraphCache *cache);
void graph_cache_reset(GraphCache *cache);

#define graph_cache_node(cache,nodeid) (&(cache)->node_buf.b[nodeid])
#define graph_cache_snode(cache,snodeid) (&(cache)->snode_buf.b[snodeid])
#define graph_cache_step(cache,stepid) (&(cache)->step_buf.b[stepid])
#define graph_cache_path(cache,pathid) (&(cache)->path_buf.b[pathid])

#define graph_cache_num_nodes(cache) ((cache)->node_buf.len)
#define graph_cache_num_snodes(cache) ((cache)->snode_buf.len)
#define graph_cache_num_steps(cache) ((cache)->step_buf.len)
#define graph_cache_num_paths(cache) ((cache)->path_buf.len)

#define graph_cache_first_node(cache,snode) \
        graph_cache_node(cache, (snode)->first_node_id)

#define graph_cache_last_node(cache,snode) \
        graph_cache_node(cache, (snode)->first_node_id + (snode)->num_nodes - 1)

#define graph_cache_path_last_step(cache,path) \
        graph_cache_step(cache, (path)->first_step + (path)->num_steps - 1)

// Returns pathid
uint32_t graph_cache_new_path(GraphCache *cache);

// Returns stepid
uint32_t graph_cache_new_step(GraphCache *cache, dBNode node);

//
// Sorting
//

int graph_cache_pathids_cmp(const void *aa, const void *bb, void *arg);

int graph_cache_steps_cmp(const GCacheStep *a, const GCacheStep *b,
                          const GraphCache *cache);

void graph_cache_steps_qsort(GraphCache *cache, GCacheStep **list, size_t n);

bool graph_cache_pathids_are_equal(GraphCache *cache,
                                   uint32_t pathid0, uint32_t pathid1);



// Get all nodes in a single step (supernode with orientation)
// Adds to the end of the node buffer (does not reset it)
void graph_cache_snode_fetch_nodes(const GraphCache *cache,
                                   const GCacheSnode *snode,
                                   Orientation orient,
                                   dBNodeBuffer *nbuf);

void graph_cache_path_fetch_nodes(const GraphCache *cache,
                                  const GCachePath *path, size_t num_steps,
                                  dBNodeBuffer *nbuf);

// Get all nodes in a path up to, but not including the given step
// Adds to the end of the node buffer (does not reset it)
void graph_cache_step_fetch_nodes(const GraphCache *cache,
                                  const GCacheStep *end_step,
                                  dBNodeBuffer *nbuf);

// Looks like 3p flank if steps don't have the same n-1 supernode
bool graph_cache_is_3p_flank(GraphCache *cache,
                             GCacheStep ** steps, size_t num_steps);

// Remove duplicate paths
size_t graph_cache_remove_dupes(GraphCache *cache,
                                GCacheStep **steps, size_t num_steps);

// Returns true if all nodes in supernode have given colour
bool graph_cache_snode_has_colour(const GraphCache *cache,
                                  const GCacheSnode *snode,
                                  size_t colour);

// Returns true if all nodes in path have given colour
bool graph_cache_step_has_colour(const GraphCache *cache,
                                 const GCacheStep *endstep,
                                 size_t colour);

// Returns NULL if not found
GCacheSnode* graph_cache_find_snode(GraphCache *cache, dBNode node);

Orientation graph_cache_get_supernode_orient(const GraphCache *cache,
                                             const GCacheSnode *snode,
                                             dBNode first_node);

#endif /* GRAPH_CACHE_H_ */
