#ifndef LINK_TREE_H_
#define LINK_TREE_H_

#include "common_buffers.h"
#include "madcrowlib/madcrow_buffer.h"

/*
 * A link tree can only represent a single sample -- they cannot be multicoloured
 * Node in the tree represent junctions in the graph. They can have up to four
 * out edges (A,C,G,T).
 */

typedef int32_t LTreeID;

typedef struct
{
  const LTreeID id, parentid; // parent is -1 if root
  LTreeID children[4];
  size_t counts[4];
  uint32_t dist; // distance to root in kmers
  uint32_t seq; // coverted to char*, sequence before this junction: dist+1 bases
  int8_t base; // choice to get here 0..3, -1 if root
} LinkJunction;

madcrow_buffer(lj_buf, LJBuffer, LinkJunction);

// Tree walking
typedef struct {
  LinkJunction *parent;
  uint8_t nxt; // 0 - 4
} LTreeWalk;

madcrow_buffer(ltree_walk_buf, LTreeWalkBuffer, LTreeWalk);

// Link Tree
typedef struct
{
  const size_t kmer_size;
  LJBuffer treebuf;
  ByteBuffer seqbuf;
  LTreeID fw_id, rv_id;
  LTreeWalkBuffer wbuf; // iterator
} LinkTree;

#define ltree_get_node(tree,id) ((id) < 0 ? NULL : (tree)->treebuf.b + (id))
#define ltree_get_fw_node(tree) ltree_get_node(tree,(tree)->fw_id)
#define ltree_get_rv_node(tree) ltree_get_node(tree,(tree)->rv_id)
#define ltree_get_seq(tree,l) ((char*)((tree)->seqbuf.b + (l)->seq))

//
// LinkTree setup
//
void ltree_alloc(LinkTree *tree, size_t kmer_size);
void ltree_dealloc(LinkTree *tree);
void ltree_reset(LinkTree *tree);

void ltree_add(LinkTree *tree,
               bool fw, size_t covg, size_t *dists,
               const char *juncs, const char *seq);

//
// Whole tree operations
//
typedef struct {
  size_t num_trees_with_links, num_links, num_link_bytes;
} LinkTreeStats;

/**
 * Get stats about this LinkTree. Does not reset LinkTreeStats.
 */
void ltree_get_stats(LinkTree *tree, LinkTreeStats *stats);

void ltree_clean(LinkTree *tree, size_t cutoff);

void ltree_write_list(LinkTree *tree, StrBuf *sbuf);

// Get number of links with ltree_get_stats() first
void ltree_write_ctp(LinkTree *tree, const char *kmer, size_t num_links,
                     StrBuf *sbuf);

void ltree_write_dot(LinkTree *tree, StrBuf *sbuf);

/**
 * @param hists should be of size distsize * covgsize
 */
void ltree_update_covg_hists(LinkTree *tree, uint64_t *hists,
                             size_t distsize, size_t covgsize);

#endif /* LINK_TREE_H_ */
