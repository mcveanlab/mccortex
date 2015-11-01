#ifndef GRAPH_CRAWLER_H_
#define GRAPH_CRAWLER_H_

#include "graph_cache.h"
#include "graph_walker.h"
#include "repeat_walker.h"

typedef struct {
  int pathid; // set to -1 if not used
  uint32_t colour;
} GCUniColPath;

typedef struct
{
  uint32_t pathid;
  uint32_t *cols, num_cols;
} GCMultiColPath;

typedef struct
{
  uint32_t num_paths;
  int *col_paths; // one per colour, col_paths[i] = -1 if colour i has no path
  GCMultiColPath *multicol_paths; // one path with multiple colours

  // used internally
  GCUniColPath *unicol_paths;
  uint32_t *col_list;
  GraphCache cache;

  // Temporary variables for walking
  GraphWalker wlk;
  RepeatWalker rptwlk;
} GraphCrawler;

static inline uint32_t graph_crawler_path_colour(const GraphCrawler *c,
                                                 uint32_t pathid,
                                                 uint32_t idx)
{
  return (c->multicol_paths[pathid].cols[idx]);
}

static inline uint32_t graph_crawler_path_ncols(const GraphCrawler *c,
                                                uint32_t pathid)
{
  return (c->multicol_paths[pathid].num_cols);
}

void graph_crawler_alloc(GraphCrawler *crawler, const dBGraph *db_graph);
void graph_crawler_dealloc(GraphCrawler *crawler);

// You don't have to reset the crawler but it wipes the cache which will
// reduce memory usage
static inline void graph_crawler_reset(GraphCrawler *crawler) {
  graph_cache_reset(&crawler->cache);
}


/**
 * @param node1 should be the first node of a unitig
 * @param node0 should be the previous node
 * @param next_base is the last base of `node1`
 * @param jmpfunc is called with each unitig traversed and if it returns true
 *                we continue crawling, otherwise we stop. If NULL assume always true
 * @param endfunc is a function called at the end of traversal
 */
void graph_crawler_fetch(GraphCrawler *crawler, dBNode node0,
                         dBNode next_nodes[4],
                         size_t take_idx, size_t num_next,
                         uint32_t *cols, size_t ncols,
                         bool (*jmpfunc)(const GraphCache *_c, const GCacheStep *_s, void *_a),
                         void (*endfunc)(const GraphCache *_c, uint32_t _pathid, void *_a),
                         void *arg);

void graph_crawler_get_path_nodes(const GraphCrawler *crawler, size_t pidx,
                                  dBNodeBuffer *nbuf);

//
// General functions on GraphCache using GraphWalker and RepeatWalker
//

/**
 * Constructs a path of unitigs (GCachePath)
 * @param wlk GraphWalker should be set to go at `node`
 * @param rptwlk RepeatWalker should be clear
 * @param jmpfunc is called with each unitig traversed and if it returns true
                  we continue crawling, otherwise we stop.
                  If NULL assume always true
 * @return pathid in GraphCache
 */
uint32_t graph_crawler_load_path(GraphCache *cache, dBNode node,
                                 GraphWalker *wlk, RepeatWalker *rptwlk,
                                 bool (*jmpfunc)(const GraphCache *_c,
                                                 const GCacheStep *_s, void *_a),
                                 void *arg);

/**
 * Constructs a path of unitigs (GCachePath)
 * @param wlk GraphWalker should be set to go at `node`
 * @param rptwlk RepeatWalker should be clear
 * @param kmer_length_limit is the max length to crawl
 * @return pathid in GraphCache
 **/
uint32_t graph_crawler_load_path_limit(GraphCache *cache, dBNode node,
                                       GraphWalker *wlk, RepeatWalker *rptwlk,
                                       size_t kmer_length_limit);

// Remove traversal marks in RepeatWalker after walking along a given path
void graph_crawler_reset_rpt_walker(RepeatWalker *rptwlk,
                                    const GraphCache *cache, uint32_t pathid);

// data[0] is number of kmers so far
// data[1] is the kmer limit
static inline bool gcrawler_load_path_limit_kmer_len(const GraphCache *cache,
                                                     const GCacheStep *step,
                                                     void *arg)
{
  const GCacheUnitig *unitig = gc_step_get_unitig(cache, step);
  size_t *data = (size_t*)arg;
  data[0] += unitig->num_nodes;
  return (data[0] < data[1]);
}

#endif /* GRAPH_CRAWLER_H_ */
