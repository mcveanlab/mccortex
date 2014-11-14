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

#define graph_crawler_path_colour(crawler,path,i) \
        ((crawler)->multicol_paths[path].cols[i])

#define graph_crawler_path_ncols(crawler,path) \
        ((crawler)->multicol_paths[path].num_cols)

void graph_crawler_alloc(GraphCrawler *crawler, const dBGraph *db_graph);
void graph_crawler_dealloc(GraphCrawler *crawler);

// You don't have to reset the crawler but it wipes the cache which will
// reduce memory usage
#define graph_crawler_reset(crawler) graph_cache_reset(&(crawler)->cache)

/**
 * @param node1 should be the first node of a supernode
 * @param node0 should be the previous node
 * @param next_base is the last base of `node1`
 * @param jmpfunc is called with each supernode traversed and if it returns true
 *                we continue crawling, otherwise we stop. If NULL assume always true
 * @param endfunc is a function called at the end of traversal
 */
void graph_crawler_fetch(GraphCrawler *crawler, dBNode node0,
                         dBNode next_nodes[4],
                         size_t take_idx, size_t num_next,
                         uint32_t *cols, size_t ncols,
                         bool (*jmpfunc)(GraphCache *_cache, GCacheStep *_step, void *_arg),
                         void (*endfunc)(GraphCache *_cache, uint32_t _pathid, void *_arg),
                         void *arg);

void graph_crawler_get_path_nodes(const GraphCrawler *crawler, size_t pidx,
                                  dBNodeBuffer *nbuf);

//
// General functions on GraphCache using GraphWalker and RepeatWalker
//

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
                                  void *arg);

// Constructs a path of supernodes (SupernodePath)
// `wlk` GraphWalker should be set to go at `node`
// `rptwlk` RepeatWalker should be clear
// `kmer_length_limit` is the max length to crawl
// returns pathid in GraphCache
uint32_t graph_crawler_load_path_limit(GraphCache *cache, dBNode node,
                                       GraphWalker *wlk, RepeatWalker *rptwlk,
                                       size_t kmer_length_limit);

// Remove traversal marks in RepeatWalker after walking along a given path
void graph_crawler_reset_rpt_walker(RepeatWalker *rptwlk,
                                    const GraphCache *cache, uint32_t pathid);

// data[0] is number of kmers so far
// data[1] is the kmer limit
static inline bool gcrawler_load_path_limit_kmer_len(GraphCache *cache,
                                                     GCacheStep *step,
                                                     void *arg)
{
  const GCacheSnode *snode = graph_cache_snode(cache, step->supernode);
  size_t *data = (size_t*)arg;
  data[0] += snode->num_nodes;
  return (data[0] < data[1]);
}

#endif /* GRAPH_CRAWLER_H_ */
