#ifndef REPEAT_WALKER_H_
#define REPEAT_WALKER_H_

#include "graph_walker.h"
#include "db_node.h"

typedef struct
{
  uint64_t *const visited, *const bloom;
  const size_t bloom_nbits, mem_bytes;
  const uint32_t mask;
  size_t nbloom_entries;
} RepeatWalker;

// GraphWalker wlk is proposing node and orient as next move
// We determine if it is safe to make the traversal without getting stuck in
// a loop/cycle in the graph
static inline bool rpt_walker_attempt_traverse(RepeatWalker *rpt,
                                               GraphWalker *wlk)
{
  uint64_t h[3], hash64;
  bool collision;

  if(!db_node_has_traversed(rpt->visited, wlk->node)) {
    db_node_set_traversed(rpt->visited, wlk->node);
    return true;
  }
  else {
    hash64 = graph_walker_hash64(wlk);

    h[0] = hash64 & rpt->mask;
    h[1] = (hash64 >>  rpt->bloom_nbits) & rpt->mask;
    h[2] = (hash64 >> (rpt->bloom_nbits*2)) & rpt->mask;

    // printf(" %zu %zu %zu / %zu\n", (size_t)h[0], (size_t)h[1], (size_t)h[2],
    //        (size_t)(1UL<<rpt->bloom_nbits));

    collision = bitset_get(rpt->bloom, h[0]) &&
                bitset_get(rpt->bloom, h[1]) &&
                bitset_get(rpt->bloom, h[2]);

    bitset_set(rpt->bloom, h[0]);
    bitset_set(rpt->bloom, h[1]);
    bitset_set(rpt->bloom, h[2]);

    rpt->nbloom_entries++;
    return !collision;
  }
}

static inline size_t rpt_walker_est_mem(size_t hash_capacity, size_t nbits)
{
  size_t visited_words = roundup_bits2words64(hash_capacity*2);
  size_t repeat_words = roundup_bits2words64(1UL<<nbits);
  return (visited_words+repeat_words) * sizeof(uint64_t);
}

static inline void rpt_walker_alloc(RepeatWalker *rpt,
                                    size_t hash_capacity, size_t nbits)
{
  ctx_assert(nbits > 0 && nbits < 32);
  size_t visited_words = roundup_bits2words64(hash_capacity*2);
  size_t repeat_words = roundup_bits2words64(1UL<<nbits);
  size_t nbytes = (visited_words + repeat_words) * sizeof(uint64_t);
  uint64_t *mem = calloc2(visited_words+repeat_words, sizeof(uint64_t));
  uint32_t mask = bitmask(nbits,uint32_t);
  RepeatWalker tmp = {.visited = mem, .bloom = mem+visited_words,
                      .bloom_nbits = nbits, .mem_bytes = nbytes, .mask = mask,
                      .nbloom_entries = 0};
  memcpy(rpt, &tmp, sizeof(RepeatWalker));
}

static inline void rpt_walker_dealloc(RepeatWalker *rpt)
{
  free(rpt->visited);
}

static inline void rpt_walker_clear(RepeatWalker *rpt)
{
  memset(rpt->visited, 0, rpt->mem_bytes);
  rpt->nbloom_entries = 0;
}

static inline void rpt_walker_fast_clear(RepeatWalker *rpt,
                                         const dBNode *nodes, size_t n)
{
  size_t i, bmem = roundup_bits2words64(1UL<<rpt->bloom_nbits)*sizeof(uint64_t);
  for(i = 0; i < n; i++) db_node_fast_clear_traversed(rpt->visited, nodes[i].key);
  if(rpt->nbloom_entries) memset(rpt->bloom, 0, bmem);
  rpt->nbloom_entries = 0;
}

static inline void rpt_walker_fast_clear_single_node(RepeatWalker *rpt,
                                                     const dBNode node)
{
  db_node_fast_clear_traversed(rpt->visited, node.key);
}

#endif /* REPEAT_WALKER_H_ */
