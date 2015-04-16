#ifndef LINK_TREE_H_
#define LINK_TREE_H_

#include "common_buffers.h"
#include "madcrowlib/madcrow_buffer.h"

typedef int32_t ltree_id;

typedef struct
{
  const int32_t id, parentid;
  int32_t children[4]; // parent is -1 if root
  uint32_t dist; // distance to root in kmers
  size_t seq; // coverted to char*, sequence before this junction: dist+1 bases
  char base; // choice to get here, -1 if root
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
  const size_t ncols, kmer_size;
  LJBuffer treebuf;
  SizeBuffer covgbuf;
  ByteBuffer seqbuf;
  int32_t fw_id, rv_id;
  int32_t fw_kstr_idx, rv_kstr_idx; // index of seqbuf of kmerstr
  LTreeWalkBuffer wbuf; // iterator
} LinkTree;

#define ltree_get_fw_kmer(tree) ((char*)((tree)->seqbuf.b + (tree)->fw_kstr_idx))
#define ltree_get_rv_kmer(tree) ((char*)((tree)->seqbuf.b + (tree)->rv_kstr_idx))
#define ltree_get_fw_node(tree) ltree_get_node(tree,(tree)->fw_id)
#define ltree_get_rv_node(tree) ltree_get_node(tree,(tree)->rv_id)
#define ltree_get_node(tree,id) ((id) < 0 ? NULL : &(tree)->treebuf.b[id])
#define ltree_get_covg(tree,id,col) ((tree)->covgbuf.b[(id)*(tree)->ncols+(col)])
#define ltree_get_seq(tree,l) ((char*)((tree)->seqbuf.b + (l)->seq))
#define ltree_node_is_leaf(l) ((l)->children[0] < 0 && (l)->children[1] < 0 && \
                               (l)->children[2] < 0 && (l)->children[3])

void ltree_remove_node(LinkTree *tree, LinkJunction *node);
size_t ltree_get_child_covg(const LinkJunction *l, const LinkTree *tree, size_t col);
// Returns true if a colour stops in the given node `l`
bool ltree_node_is_col_leaf(const LinkJunction *l, const LinkTree *tree);
bool ltree_node_has_covg(const LinkJunction *l, const LinkTree *tree);

#define ltree_link_col_covg(tree,lj,col) (ltree_get_covg(tree,(lj)->id,col) - ltree_get_child_covg(l,tree,col))

//
// LinkTree setup
//
void ltree_alloc(LinkTree *tree, size_t ncols, size_t kmer_size);
void ltree_dealloc(LinkTree *tree);
void ltree_reset(LinkTree *tree);

void ltree_add(LinkTree *tree,
               bool fw, size_t *covgs, size_t *dists,
               const char *juncs, const char *seq);

//
// Tree walking
//
LinkJunction* ltree_visit_nodes_sub(LinkTree *tree, LinkJunction *root,
                                    bool (*func)(LinkJunction *_lj,
                                                 LinkTree *_tree,
                                                 uint32_t _depth,
                                                 void *_ptr),
                                    void *ptr);

void ltree_visit_nodes(LinkTree *tree,
                       bool (*func)(LinkJunction *_lj,
                                    LinkTree *_tree,
                                    uint32_t _depth,
                                    void *_ptr),
                       void *ptr);

//
// Whole tree operations
//
typedef struct {
  size_t num_trees_with_links, num_links, num_link_bytes;
} LinkTreeStats;

// Returns true if this tree contains links
bool ltree_count_col_leaves(LinkTree *tree, LinkTreeStats *stats);

void ltree_clean(LinkTree *tree, size_t *cutoffs);

void ltree_write_list(LinkTree *tree, StrBuf *sbuf);

// Get number of links with ltree_count_col_leaves() first
void ltree_write_ctp(LinkTree *tree, const char *kmer, size_t num_links,
                     StrBuf *sbuf);

void ltree_write_dot(LinkTree *tree, bool fw, const char *kmer, StrBuf *sbuf);

void ltree_update_covg_hists(LinkTree *tree, uint64_t *hists,
                             size_t distlen, size_t covglen);

#endif /* LINK_TREE_H_ */
