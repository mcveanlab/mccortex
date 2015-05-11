#include "global.h"
#include "graph_crawler.h"
#include "binary_seq.h" // binary_seq_unpack_byte()
#include "sort_r/sort_r.h" // binary_seq_unpack_byte()

static inline void walk_supernode_end(const GraphCache *cache,
                                      const GCacheSnode *snode,
                                      Orientation snorient,
                                      GraphWalker *wlk)
{
  // Only need to traverse the first and last nodes of a supernode
  const dBNode *first = graph_cache_first_node(cache, snode);
  const dBNode *last = graph_cache_last_node(cache, snode);
  dBNode lastnode;
  size_t num_nodes = last - first + 1;

  if(num_nodes > 1) {
    lastnode = (snorient == FORWARD ? *last : db_node_reverse(*first));
    graph_walker_jump_along_snode(wlk, lastnode, num_nodes);
  }
}

// Constructs a path of supernodes (SupernodePath)
// `wlk` GraphWalker should be set to go at `node`
// `rptwlk` RepeatWalker should be clear
// `jmpfunc` is called with each supernode traversed and if it returns true
//           we continue crawling, otherwise we stop. If NULL assume always true
// returns pathid in GraphCache
uint32_t graph_crawler_load_path(GraphCache *cache, dBNode node,
                                 GraphWalker *wlk, RepeatWalker *rptwlk,
                                 bool (*jmpfunc)(GraphCache *_cache,
                                                 GCacheStep *_step, void *_arg),
                                 void *arg)
{
  size_t i;
  uint32_t stepid, pathid = graph_cache_new_path(cache);

  ctx_assert(db_nodes_are_equal(wlk->node, node));

  for(i = 0; ; i++)
  {
    stepid = graph_cache_new_step(cache, node);

    GCacheStep *step = graph_cache_step(cache, stepid);
    GCacheSnode *snode = graph_cache_snode(cache, step->supernode);

    // Traverse to the end of the supernode
    walk_supernode_end(cache, snode, step->orient, wlk);

    if(jmpfunc != NULL && !jmpfunc(cache, step, arg)) break;

    // Find next node
    uint8_t num_edges;
    const dBNode *next_nodes;
    Nucleotide next_bases[4];

    if(step->orient == FORWARD) {
      num_edges = snode->num_next;
      next_nodes = snode->next_nodes;
      binary_seq_unpack_byte(next_bases, snode->next_bases);
    }
    else {
      num_edges = snode->num_prev;
      next_nodes = snode->prev_nodes;
      binary_seq_unpack_byte(next_bases, snode->prev_bases);
    }

    // Traverse to next supernode
    if(!graph_walker_next_nodes(wlk, num_edges, next_nodes, next_bases) ||
       !rpt_walker_attempt_traverse(rptwlk, wlk)) break;

    node = wlk->node;
  }

  return pathid;
}

// Constructs a path of supernodes (SupernodePath)
// `wlk` GraphWalker should be set to go at `node`
// `rptwlk` RepeatWalker should be clear
// returns pathid in GraphCache
uint32_t graph_crawler_load_path_limit(GraphCache *cache, dBNode node,
                                       GraphWalker *wlk, RepeatWalker *rptwlk,
                                       size_t kmer_length_limit)
{
  size_t path_limit[2] = {0, kmer_length_limit};

  return graph_crawler_load_path(cache, node, wlk, rptwlk,
                                 gcrawler_load_path_limit_kmer_len, path_limit);
}


void graph_crawler_reset_rpt_walker(RepeatWalker *rptwlk,
                                    const GraphCache *cache, uint32_t pathid)
{
  rpt_walker_fast_clear(rptwlk, NULL, 0);

  const GCachePath *path = graph_cache_path(cache, pathid);
  const GCacheStep *step = graph_cache_step(cache, path->first_step), *endstep;
  const GCacheSnode *snode;
  const dBNode *node0, *node1;

  // Loop over supernodes in the path
  for(endstep = step + path->num_steps; step < endstep; step++)
  {
    // We don't care about orientation here
    snode = graph_cache_snode(cache, step->supernode);
    node0 = graph_cache_first_node(cache, snode);
    node1 = graph_cache_last_node(cache, snode);
    rpt_walker_fast_clear_single_node(rptwlk, *node0);
    rpt_walker_fast_clear_single_node(rptwlk, *node1);
  }
}


void graph_crawler_alloc(GraphCrawler *crawler, const dBGraph *db_graph)
{
  ctx_assert(db_graph->node_in_cols != NULL);

  size_t ncols = db_graph->num_of_cols;

  int *col_paths = ctx_calloc(ncols, sizeof(int));
  GCMultiColPath *multicol_paths = ctx_calloc(ncols, sizeof(GCMultiColPath));
  GCUniColPath *unicol_paths = ctx_calloc(ncols, sizeof(GCUniColPath));
  uint32_t *col_list = ctx_calloc(ncols, sizeof(uint32_t));

  GraphCrawler tmp = {.num_paths = 0,
                      .col_paths = col_paths,
                      .multicol_paths = multicol_paths,
                      .unicol_paths = unicol_paths,
                      .col_list = col_list};

  memcpy(crawler, &tmp, sizeof(GraphCrawler));

  graph_cache_alloc(&crawler->cache, db_graph);
  graph_walker_alloc(&crawler->wlk, db_graph);
  rpt_walker_alloc(&crawler->rptwlk, db_graph->ht.capacity, 22); // 4MB
}

void graph_crawler_dealloc(GraphCrawler *crawler)
{
  ctx_free(crawler->col_paths);
  ctx_free(crawler->multicol_paths);
  ctx_free(crawler->unicol_paths);
  ctx_free(crawler->col_list);
  graph_cache_dealloc(&crawler->cache);
  graph_walker_dealloc(&crawler->wlk);
  rpt_walker_dealloc(&crawler->rptwlk);
  memset(crawler, 0, sizeof(GraphCrawler)); // reset
}

static inline int unicol_path_cmp(const void *aa, const void *bb, void *arg)
{
  const GCMultiColPath *a = (const GCMultiColPath *)aa;
  const GCMultiColPath *b = (const GCMultiColPath *)bb;
  GraphCache *cache = (GraphCache*)arg;
  return graph_cache_pathids_cmp(&a->pathid, &b->pathid, cache);
}

// `node1` should be the first node of a supernode
// `node0` should be the previous node
// `next_base` is the last base of `node1`
// `jmpfunc` is called with each supernode traversed and if it returns true
//           we continue crawling, otherwise we stop
// `endfunc` is a function called at the end of traversal
void graph_crawler_fetch(GraphCrawler *crawler, dBNode node0,
                         dBNode next_nodes[4],
                         size_t take_idx, size_t num_next,
                         uint32_t *cols, size_t ncols,
                         bool (*jmpfunc)(GraphCache *_cache, GCacheStep *_step, void *_arg),
                         void (*endfunc)(GraphCache *_cache, uint32_t _pathid, void *_arg),
                         void *arg)
{
  const dBGraph *db_graph = crawler->cache.db_graph;
  GraphCache *cache = &crawler->cache;
  GraphWalker *wlk = &crawler->wlk;
  RepeatWalker *rptwlk = &crawler->rptwlk;
  GCUniColPath *unipaths = crawler->unicol_paths;

  ctx_assert(take_idx < num_next);
  ctx_assert(!db_nodes_are_equal(node0, next_nodes[take_idx]));

  // Fetch all paths in all colours
  dBNode node1 = next_nodes[take_idx];
  bool is_fork;
  size_t i, c, col, nedges_cols, num_unicol_paths = 0;
  int pathid;

  for(c = 0; c < ncols; c++)
  {
    col = (cols != NULL ? cols[c] : c);

    if(db_node_has_col(db_graph, node0.key, col) &&
       db_node_has_col(db_graph, node1.key, col))
    {
      // Determine if this fork is a fork in the current colour
      for(nedges_cols = 0, i = 0; i < num_next && nedges_cols <= 1; i++)
        nedges_cols += db_node_has_col(db_graph, next_nodes[i].key, col);

      is_fork = (nedges_cols > 1);

      graph_walker_setup(wlk, true, col, col, db_graph);
      graph_walker_start(wlk, node0);
      graph_walker_force(wlk, node1, is_fork);

      pathid = graph_crawler_load_path(cache, node1, wlk, rptwlk, jmpfunc, arg);

      if(endfunc != NULL) endfunc(cache, pathid, arg);

      graph_walker_finish(wlk);
      graph_crawler_reset_rpt_walker(rptwlk, cache, pathid);

      unipaths[num_unicol_paths++] = (GCUniColPath){.colour = col,
                                                    .pathid = pathid};
    }
    else
      pathid = -1;

    crawler->col_paths[col] = pathid;
  }

  if(num_unicol_paths == 0) {
    crawler->num_paths = 0;
  }
  else {
  // sort unicol paths to group duplicate paths
    sort_r(unipaths, num_unicol_paths, sizeof(GCUniColPath), unicol_path_cmp, cache);

    // Create Multicol paths by merging adjacent identical unicol paths
    uint32_t *col_list = crawler->col_list;
    GCMultiColPath *multicol_paths = crawler->multicol_paths;

    for(i = 0; i < num_unicol_paths; i++) col_list[i] = unipaths[i].colour;

    size_t num_multipaths = 1;
    pathid = unipaths[0].pathid;
    multicol_paths[0] = (GCMultiColPath){.pathid = pathid,
                                         .cols = col_list,
                                         .num_cols = 1};

    for(i = 1; i < num_unicol_paths; i++)
    {
      if(graph_cache_pathids_are_equal(cache, unipaths[i].pathid, pathid)) {
        // Path matches existing
        multicol_paths[num_multipaths-1].num_cols++;
      }
      else {
        // New path
        col_list += multicol_paths[num_multipaths-1].num_cols;
        pathid = unipaths[i].pathid;
        multicol_paths[num_multipaths] = (GCMultiColPath){.pathid = pathid,
                                                          .cols = col_list,
                                                          .num_cols = 1};
        num_multipaths++;
      }
    }

    ctx_assert(num_multipaths <= num_unicol_paths);
    ctx_assert(num_multipaths <= db_graph->num_of_cols);

    crawler->num_paths = num_multipaths;
  }
}

void graph_crawler_get_path_nodes(const GraphCrawler *crawler, size_t pidx,
                                  dBNodeBuffer *nbuf)
{
  uint32_t pathid = crawler->multicol_paths[pidx].pathid;
  const GraphCache *cache = &crawler->cache;
  const GCachePath *path = graph_cache_path(cache, pathid);
  graph_cache_path_fetch_nodes(cache, path, path->num_steps, nbuf);
}
