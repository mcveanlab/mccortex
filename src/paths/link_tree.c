#include "global.h"
#include "link_tree.h"
#include "dna.h"

void ltree_remove_node(LinkTree *tree, LinkJunction *node)
{
  if(node->parentid > 0)
    ltree_get_node(tree, node->parentid)->children[(size_t)node->base] = -1;
  else if(tree->fw_id == node->id)
    tree->fw_id = -1;
  else if(tree->rv_id == node->id)
    tree->rv_id = -1;
  else
    die("Cannot find root node in tree [%i; fw:%i; rv:%i]",
        node->id, tree->fw_id, tree->rv_id);
}

size_t ltree_get_child_covg(const LinkJunction *l, const LinkTree *tree, size_t col)
{
  size_t i, sum;
  for(i = sum = 0; i < 4; i++)
    if(l->children[i] >= 0)
      sum += ltree_get_covg(tree,l->children[i],col);
  return sum;
}

bool ltree_node_is_col_leaf(const LinkJunction *l, const LinkTree *tree)
{
  size_t col;

  for(col = 0; col < tree->ncols; col++) {
    if(ltree_get_covg(tree,l->id,col) > 0 && ltree_get_child_covg(l,tree,col) == 0)
      return true;
  }

  return false;
}

bool ltree_node_has_covg(const LinkJunction *l, const LinkTree *tree)
{
  size_t col, total = 0;
  for(col = 0; col < tree->ncols; col++)
    total += ltree_get_covg(tree, l->id, col);
  return (total > 0);
}

static inline bool _count_col_leaves(LinkJunction *l, LinkTree *tree,
                                     uint32_t depth, void *ptr)
{
  size_t njuncs = depth+1;
  LinkTreeStats *stats = (LinkTreeStats*)ptr;
  if(ltree_node_is_col_leaf(l,tree)) {
    stats->num_links++;
    stats->num_link_bytes += (njuncs+3)/4;
  }
  return true;
}

// Returns true if this tree contains links
bool ltree_count_col_leaves(LinkTree *tree, LinkTreeStats *stats)
{
  size_t init_num_links = stats->num_links;
  ltree_visit_nodes(tree, _count_col_leaves, stats);
  bool tree_has_links = (init_num_links < stats->num_links);
  stats->num_trees_with_links += tree_has_links;
  return tree_has_links;
}

//
// LinkTree setup
//

void ltree_alloc(LinkTree *tree, size_t ncols, size_t kmer_size)
{
  LinkTree tmp = {.ncols = ncols, .kmer_size = kmer_size,
                  .fw_id = -1, .rv_id = -1};
  memcpy(tree, &tmp, sizeof(tmp));
  lj_buf_alloc(&tree->treebuf, 128);
  size_buf_alloc(&tree->covgbuf, 128*ncols);
  ltree_walk_buf_alloc(&tree->wbuf, 128);
}

void ltree_dealloc(LinkTree *tree)
{
  lj_buf_dealloc(&tree->treebuf);
  size_buf_dealloc(&tree->covgbuf);
  ltree_walk_buf_dealloc(&tree->wbuf);
  memset(tree, 0, sizeof(*tree));
  tree->fw_id = tree->rv_id = -1;
}

void ltree_reset(LinkTree *tree)
{
  lj_buf_reset(&tree->treebuf);
  size_buf_reset(&tree->covgbuf);
  ltree_walk_buf_reset(&tree->wbuf);
  tree->fw_id = tree->rv_id = -1;
}

static int32_t ltree_init_node(LinkTree *tree, int32_t parent, uint32_t dist,
                               const char *seq, size_t seqlen,
                               char base)
{
  size_t col;
  LinkJunction tmp = {.id = tree->treebuf.len,
                      .parentid = parent, .children = {-1,-1,-1,-1},
                      .dist = dist, .seq = tree->seqbuf.len,
                      .base = base};
  lj_buf_add(&tree->treebuf, tmp);
  byte_buf_append(&tree->seqbuf, (const uint8_t*)seq, seqlen);
  byte_buf_add(&tree->seqbuf, '\0');
  for(col = 0; col < tree->ncols; col++) size_buf_add(&tree->covgbuf, 0);
  ctx_assert((tmp.id+1) * tree->ncols == tree->covgbuf.len);
  return tmp.id;
}

void ltree_add(LinkTree *tree,
               bool fw, size_t *covgs, size_t *dists,
               const char *juncs, const char *seq)
{
  size_t i, col, njuncs = strlen(juncs);
  int32_t rootid = fw ? tree->fw_id : tree->rv_id;

  // Sanity check: seq should include the kmer
  ctx_assert(strlen(seq) == tree->kmer_size + dists[njuncs-1] + 1);

  if(rootid < 0) {
    rootid = ltree_init_node(tree, -1, dists[0],
                             seq, tree->kmer_size + dists[0],
                             seq[tree->kmer_size+dists[0]]);
    if(fw) tree->fw_id = rootid;
    else   tree->rv_id = rootid;

    for(col = 0; col < tree->ncols; col++)
      ltree_get_covg(tree,rootid,col) += covgs[col];
  }

  int32_t parentid = rootid, nodeid;

  for(i = 1; i < njuncs; i++, parentid = nodeid)
  {
    char base = seq[tree->kmer_size+dists[i]];
    size_t baseidx = dna_char_to_nuc(base);

    nodeid = ltree_get_node(tree, parentid)->children[baseidx];
    if(nodeid < 0) {
      nodeid = ltree_init_node(tree, parentid, dists[i],
                               seq + tree->kmer_size + dists[i],
                               dists[i] - dists[i-1] - 1,
                               base);
      ltree_get_node(tree, parentid)->children[baseidx] = nodeid;
    }

    for(col = 0; col < tree->ncols; col++)
      ltree_get_covg(tree,nodeid,col) += covgs[col];
  }
}


//
// Tree Walking
//

LinkJunction* ltree_visit_nodes_sub(LinkTree *tree, LinkJunction *root,
                                    bool (*func)(LinkJunction *_lj,
                                                 LinkTree *_tree,
                                                 uint32_t _depth,
                                                 void *_ptr),
                                    void *ptr)
{
  if(!root || !func(root, tree, 0, ptr)) return NULL;

  // wbuf.b[wbuf.len-1] is the next junction to take
  LTreeWalkBuffer *wbuf = &tree->wbuf;
  ltree_walk_buf_reset(wbuf);

  // take the 0th junction next time
  ltree_walk_buf_add(wbuf, (LTreeWalk){.parent = root, .nxt = 0});

  while(wbuf->len > 0)
  {
    // Attempt to take a junction
    LTreeWalk *walk = &wbuf->data[wbuf->len-1];
    while(walk->nxt < 4 && walk->parent->children[walk->nxt] < 0) walk->nxt++;

    if(walk->nxt == 4) { ltree_walk_buf_pop(wbuf); }
    else {
      LinkJunction *node = ltree_get_node(tree, walk->parent->children[walk->nxt]);

      if(func(node, tree, wbuf->len, ptr))
        ltree_walk_buf_add(wbuf, (LTreeWalk){.parent = node, .nxt = 0});
    }
  }

  return root;
}

void ltree_visit_nodes(LinkTree *tree,
                       bool (*func)(LinkJunction *_lj,
                                    LinkTree *_tree,
                                    uint32_t _depth,
                                    void *_ptr),
                       void *ptr)
{
  LinkJunction *root;
  if((root = ltree_get_fw_node(tree)) != NULL)
    ltree_visit_nodes_sub(tree, root, func, ptr);
  if((root = ltree_get_rv_node(tree)) != NULL)
    ltree_visit_nodes_sub(tree, root, func, ptr);
}

//
// Whole tree operations
//

// Returns true if we should keep traversing
static inline bool _ltree_clean_node(LinkJunction *l, LinkTree *tree,
                                     uint32_t depth, void *ptr)
{
  (void)depth;
  size_t col;
  size_t *cutoffs = (size_t*)ptr;
  bool keep_node = false;
  for(col = 0; col < tree->ncols; col++) {
    if(ltree_get_covg(tree, l->id, col) >= cutoffs[col]) keep_node = true;
    else ltree_get_covg(tree, l->id, col) = 0;
  }
  if(!keep_node) {
    // No colours have coverage
    ltree_remove_node(tree, ltree_get_node(tree,l->id));
    return false;
  }
  return true;
}

// List coverage of nodes in colour 0
static inline bool _ltree_write_list_node(LinkJunction *l, LinkTree *tree,
                                          uint32_t depth, void *ptr)
{
  (void)depth;
  size_t col;
  StrBuf *sbuf = (StrBuf*)ptr;
  strbuf_append_ulong(sbuf, tree->kmer_size + l->dist + 1);
  for(col = 0; col < tree->ncols; col++) {
    strbuf_append_char(sbuf, ',');
    strbuf_append_ulong(sbuf, ltree_get_covg(tree, l->id, col));
  }
  strbuf_append_char(sbuf, '\n');
  return true; // visit all nodes
}

static inline bool _ltree_write_dot_node(LinkJunction *l, LinkTree *tree,
                                         uint32_t depth, void *ptr)
{
  (void)depth;
  StrBuf* sbuf = (StrBuf*)ptr;
  size_t col;
  strbuf_sprintf(sbuf, "  node%i [label=\"", l->id);
  strbuf_sprintf(sbuf, "%s %c (dist: %u count: %zu",
                       ltree_get_seq(tree, l), l->base,
                       l->dist, ltree_get_covg(tree, l->id, 0));
  for(col = 1; col < tree->ncols; col++)
    strbuf_sprintf(sbuf, ",%zu", ltree_get_covg(tree, l->id, col));
  strbuf_append_str(sbuf, ")\"]\n");
  return true; // print every node
}

static inline bool _ltree_write_dot_edges(LinkJunction *l, LinkTree *tree,
                                          uint32_t depth, void *ptr)
{
  (void)tree; (void)depth;
  StrBuf* sbuf = (StrBuf*)ptr;
  size_t i;
  for(i = 0; i < 4; i++) {
    if(l->children[i] >= 0) {
      strbuf_sprintf(sbuf, "  node%i -> node%i\n", l->id, l->children[i]);
    }
  }
  return true;
}

static inline bool _ltree_write_ctp_node(LinkJunction *l, LinkTree *tree,
                                         uint32_t depth, void *ptr)
{
  (void)depth;

  // Check if a colour stops in this node
  if(!ltree_node_is_col_leaf(l,tree)) return true;

  StrBuf *sbuf = (StrBuf*)ptr;
  LinkJunction *lj;
  size_t i, col, njuncs = 0;

  // Count tree depth
  for(lj = l; lj; lj = ltree_get_node(tree,lj->parentid)) { njuncs++; }

  LinkJunction *nodes[njuncs];
  char juncs[njuncs];

  for(i = njuncs, lj = l; lj; lj = ltree_get_node(tree,lj->parentid)) {
    i--;
    nodes[i] = lj;
    juncs[i] = lj->base;
  }

  // [FR] [num_kmers] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] [seq=...] [juncpos=...]
  char dir = (nodes[0]->id == tree->fw_id ? 'F' : 'R');

  strbuf_append_char(sbuf, dir); // [FR]
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, l->dist+2); // [num_kmers]
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, njuncs); // [num_juncs]

  // [counts0,counts1,...]
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, ltree_link_col_covg(tree,l, 0));
  for(col = 1; col < tree->ncols; col++) {
    strbuf_append_char(sbuf, ',');
    strbuf_append_ulong(sbuf, ltree_link_col_covg(tree,l, col));
  }

  strbuf_append_char(sbuf, ' ');
  strbuf_append_strn(sbuf, juncs, njuncs);

  strbuf_append_str(sbuf, " seq=");
  for(i = 0; i < njuncs; i++)
    strbuf_append_str(sbuf, ltree_get_seq(tree, nodes[i]));
  strbuf_append_char(sbuf, l->base);

  strbuf_append_str(sbuf, " juncpos=");
  strbuf_append_ulong(sbuf, nodes[0]->dist);
  for(i = 1; i < njuncs; i++) {
    strbuf_append_char(sbuf, ',');
    strbuf_append_ulong(sbuf, nodes[i]->dist);
  }
  strbuf_append_char(sbuf, '\n');

  return true; // visit all nodes
}

void ltree_clean(LinkTree *tree, size_t *cutoffs)
{
  ltree_visit_nodes(tree, _ltree_clean_node, cutoffs);
}

void ltree_write_list(LinkTree *tree, StrBuf *sbuf)
{
  ltree_visit_nodes(tree, _ltree_write_list_node, sbuf);
}

// Get number of links with ltree_count_col_leaves() first
void ltree_write_ctp(LinkTree *tree, const char *kmer, size_t num_links,
                     StrBuf *sbuf)
{
  strbuf_append_str(sbuf, kmer);
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, num_links);
  strbuf_append_char(sbuf, '\n');
  ltree_visit_nodes(tree, _ltree_write_ctp_node, sbuf);
}

void ltree_write_dot(LinkTree *tree, bool fw, const char *kmer, StrBuf *sbuf)
{
  LinkJunction *root = fw ? ltree_get_fw_node(tree) : ltree_get_rv_node(tree);
  strbuf_append_str(sbuf, "digraph G {\n");
  strbuf_append_str(sbuf, "  node [shape=none fontname=\"Courier New\" fontsize=9]\n");
  strbuf_append_str(sbuf, "  kmer [label=\"");
  strbuf_append_str(sbuf, kmer);
  strbuf_append_str(sbuf, fw ? " F" : " R");
  strbuf_append_str(sbuf, "\"]\n");
  ltree_visit_nodes_sub(tree, root, _ltree_write_dot_node, sbuf);
  strbuf_append_str(sbuf, "  kmer -> node");
  strbuf_append_ulong(sbuf, root->id);
  strbuf_append_str(sbuf, "\n");
  ltree_visit_nodes_sub(tree, root, _ltree_write_dot_edges, sbuf);
  strbuf_append_str(sbuf, "}\n");
}

typedef struct {
  uint64_t *hists; // cast to [col][dist][covg]
  size_t distlen, covglen;
} LTreeHists;

static inline bool _ltree_update_covg_hists(LinkJunction *l, LinkTree *tree,
                                            uint32_t depth, void *ptr)
{
  (void)depth;
  size_t col, covg;

  LTreeHists *d = (LTreeHists*)ptr;
  uint64_t (*hists)[d->distlen][d->covglen]
    = (uint64_t (*)[d->distlen][d->covglen])d->hists;

  if(l->dist >= d->distlen) return false;

  for(col = 0; col < tree->ncols; col++) {
    covg = ltree_get_covg(tree,l->id,col);
    hists[col][l->dist][MIN2(covg,d->covglen)]++;
  }
  return true;
}

void ltree_update_covg_hists(LinkTree *tree, uint64_t *hists,
                             size_t distlen, size_t covglen)
{
  ctx_assert(hists && distlen > 0 && covglen > 0);
  LTreeHists data = {.hists = hists, .distlen = distlen, .covglen = covglen};
  ltree_visit_nodes(tree, _ltree_update_covg_hists, &data);
}
