#include "global.h"
#include "link_tree.h"
#include "dna.h"

//
// LinkTree setup
//

void ltree_alloc(LinkTree *tree, size_t kmer_size)
{
  LinkTree tmp = {.kmer_size = kmer_size, .fw_id = -1, .rv_id = -1};
  memcpy(tree, &tmp, sizeof(tmp));
  lj_buf_alloc(&tree->treebuf, 128);
  ltree_walk_buf_alloc(&tree->wbuf, 128);
}

void ltree_dealloc(LinkTree *tree)
{
  lj_buf_dealloc(&tree->treebuf);
  ltree_walk_buf_dealloc(&tree->wbuf);
  memset(tree, 0, sizeof(*tree));
  tree->fw_id = tree->rv_id = -1;
}

void ltree_reset(LinkTree *tree)
{
  lj_buf_reset(&tree->treebuf);
  byte_buf_reset(&tree->seqbuf);
  ltree_walk_buf_reset(&tree->wbuf);
  tree->fw_id = tree->rv_id = -1;
}

static LTreeID ltree_init_node(LinkTree *tree, LTreeID parent, uint32_t dist,
                               const char *seq, size_t seqlen,
                               int8_t base)
{
  LinkJunction tmp = {.id = tree->treebuf.len,
                      .parentid = parent,
                      .children = {-1,-1,-1,-1},
                      .counts = {0,0,0,0},
                      .dist = dist, .seq = tree->seqbuf.len,
                      .base = base};
  lj_buf_add(&tree->treebuf, tmp);
  byte_buf_push(&tree->seqbuf, (const uint8_t*)seq, seqlen);
  byte_buf_add(&tree->seqbuf, '\0');
  return tmp.id;
}

/**
 * Add a link to a link tree
 * @param fw    true if link starts from the kmer in the forward orientation
 * @param covg  coverage in sample (count)
 * @param dists distance from start of each junction choice (in kmers)
 * @param juncs junction choices made by this link
 * @param seq   entire sequence of the path through the graph described by link
 */
void ltree_add(LinkTree *tree,
               bool fw, size_t covg, size_t *dists,
               const char *juncs, const char *seq)
{
  size_t i, njuncs = strlen(juncs);

  // Sanity check: seq should include the kmer
  ctx_assert(strlen(seq) == tree->kmer_size + dists[njuncs-1] + 1);

  LTreeID rootid = fw ? tree->fw_id : tree->rv_id;
  LTreeID parentid = -1, nodeid = -1;
  uint32_t nuc, prev_nuc;

  for(i = 0; i < njuncs; i++, parentid = nodeid, prev_nuc = nuc)
  {
    char base = seq[tree->kmer_size+dists[i]];
    nuc = dna_char_to_nuc(base);

    nodeid = i > 0 ? ltree_get_node(tree, parentid)->children[prev_nuc] : rootid;

    if(nodeid < 0) {
      // node doesn't exist
      size_t seq_offset = i > 0 ? tree->kmer_size + dists[i-1] + 1 : 0;
      size_t dist       = i > 0 ? dists[i] - dists[i-1] - 1 : tree->kmer_size + dists[0];

      nodeid = ltree_init_node(tree, parentid, dists[i],
                               seq + seq_offset, dist,
                               i > 0 ? prev_nuc : -1);

      if(i > 0) ltree_get_node(tree, parentid)->children[prev_nuc] = nodeid;
      else if(fw) { tree->fw_id = nodeid; }
      else   { tree->rv_id = nodeid; }
    }

    // printf("%zu) nodeid: %i parentid: %i base: %i\n", i, nodeid, parentid, (int)nuc);

    ltree_get_node(tree, nodeid)->counts[nuc] += covg;
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
    LTreeWalk *walk = &wbuf->b[wbuf->len-1];
    while(walk->nxt < 4 && walk->parent->children[walk->nxt] < 0) walk->nxt++;

    if(walk->nxt == 4) { ltree_walk_buf_pop(wbuf, NULL, 1); }
    else {
      LinkJunction *node = ltree_get_node(tree, walk->parent->children[walk->nxt]);
      walk->nxt++;

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

LinkJunction* ltree_visit_links_sub(LinkTree *tree, LinkJunction *root,
                                    bool (*func)(LinkJunction *_lj, uint8_t _b,
                                                 LinkTree *_tree,
                                                 uint32_t _depth,
                                                 void *_ptr),
                                    void *ptr)
{
  if(!root) return NULL;

  // wbuf.b[wbuf.len-1] is the next junction to take
  LTreeWalkBuffer *wbuf = &tree->wbuf;
  ltree_walk_buf_reset(wbuf);

  // take the 0th junction next time
  ltree_walk_buf_add(wbuf, (LTreeWalk){.parent = root, .nxt = 0});

  while(wbuf->len > 0)
  {
    // Attempt to take a junction
    LTreeWalk *walk = &wbuf->b[wbuf->len-1];
    while(walk->nxt < 4 && walk->parent->counts[walk->nxt] == 0) walk->nxt++;
    uint8_t base = walk->nxt++; // Last base of (sub)-link

    if(base == 4) { ltree_walk_buf_pop(wbuf, NULL, 1); }
    else if(func(walk->parent, base, tree, wbuf->len-1, ptr) &&
            walk->parent->children[base] >= 0)
    {
      LinkJunction *child = ltree_get_node(tree, walk->parent->children[base]);
      ltree_walk_buf_add(wbuf, (LTreeWalk){.parent = child, .nxt = 0});
    }
  }

  return root;
}

/**
 * Iterate over all links stored in a link tree
 */
void ltree_visit_links(LinkTree *tree,
                       bool (*func)(LinkJunction *_lj, uint8_t _b,
                                    LinkTree *_tree,
                                    uint32_t _depth,
                                    void *_ptr),
                       void *ptr)
{
  LinkJunction *root;
  if((root = ltree_get_fw_node(tree)) != NULL)
    ltree_visit_links_sub(tree, root, func, ptr);
  if((root = ltree_get_rv_node(tree)) != NULL)
    ltree_visit_links_sub(tree, root, func, ptr);
}

//
// Whole tree operations
//

/* Get stats */

static inline bool _link_get_stats(LinkJunction *l, uint8_t base,
                                   LinkTree *tree, uint32_t depth, void *ptr)
{
  (void)tree;
  size_t njuncs = depth+1;
  LinkTreeStats *stats = (LinkTreeStats*)ptr;
  if(l->children[base] < 0) {
    // This link is not a substring of a longer link
    stats->num_links++;
    stats->num_link_bytes += (njuncs+3)/4;
  }
  return true; // visit all links
}

/**
 * Get stats about this LinkTree. Does not reset LinkTreeStats.
 */
void ltree_get_stats(LinkTree *tree, LinkTreeStats *stats)
{
  size_t init_num_links = stats->num_links;
  ltree_visit_links(tree, _link_get_stats, stats);
  stats->num_trees_with_links += (stats->num_links > init_num_links);
}

/* Clean links */

// Returns true if we should keep traversing
static inline bool _ltree_clean_link(LinkJunction *l, uint8_t base,
                                     LinkTree *tree,
                                     uint32_t depth, void *ptr)
{
  (void)tree; (void)depth;
  size_t cutoff = *(size_t*)ptr;
  if(l->counts[base] < cutoff) {
    l->counts[base] = 0;
    l->children[base] = -1;
    return false;
  }
  return true;
}

void ltree_clean(LinkTree *tree, size_t cutoff)
{
  ltree_visit_links(tree, _ltree_clean_link, &cutoff);
}

/* Write CSV file of link coverage */

// List coverage of links in colour 0
static inline bool _ltree_write_list_link(LinkJunction *l, uint8_t base,
                                          LinkTree *tree,
                                          uint32_t depth, void *ptr)
{
  (void)depth;
  StrBuf *sbuf = (StrBuf*)ptr;
  strbuf_append_ulong(sbuf, tree->kmer_size + l->dist + 1);
  strbuf_append_char(sbuf, ',');
  strbuf_append_ulong(sbuf, l->counts[base]);
  strbuf_append_char(sbuf, '\n');
  return true; // visit all nodes
}

void ltree_write_list(LinkTree *tree, StrBuf *sbuf)
{
  ltree_visit_links(tree, _ltree_write_list_link, sbuf);
}

/* Write .ctp file */

static inline bool _ltree_write_ctp_link(LinkJunction *l, uint8_t base,
                                         LinkTree *tree,
                                         uint32_t depth, void *ptr)
{
  // Check if this link is a prefix of a longer link
  if(l->children[base] >= 0) return true;

  StrBuf *sbuf = (StrBuf*)ptr;
  LinkJunction *lj;
  size_t i, njuncs = depth+1;

  // root + junction choices
  LinkJunction *nodes[njuncs];
  char juncs[njuncs];

  for(i = njuncs-1, lj = l; lj; lj = ltree_get_node(tree,lj->parentid), i--) {
    nodes[i] = lj;
  }

  for(i = 0; i+1 < njuncs; i++) juncs[i] = "ACGT"[nodes[i+1]->base];
  juncs[njuncs-1] = "ACGT"[base];

  // [FR] [num_kmers] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] [seq=...] [juncpos=...]
  bool fw = (nodes[0]->id == tree->fw_id);

  strbuf_append_char(sbuf, fw ? 'F' : 'R'); // [FR]
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, l->dist+2); // [num_kmers]
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, njuncs); // [num_juncs]

  // [counts0,counts1,...]
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, l->counts[base]);

  // [juncs]
  strbuf_append_char(sbuf, ' ');
  strbuf_append_strn(sbuf, juncs, njuncs);

  strbuf_append_str(sbuf, " seq=");
  for(i = 0; i < njuncs; i++) {
    strbuf_append_str(sbuf, ltree_get_seq(tree, nodes[i]));
    strbuf_append_char(sbuf, juncs[i]);
  }

  strbuf_append_str(sbuf, " juncpos=");
  strbuf_append_ulong(sbuf, nodes[0]->dist);
  for(i = 1; i < njuncs; i++) {
    strbuf_append_char(sbuf, ',');
    strbuf_append_ulong(sbuf, nodes[i]->dist);
  }
  strbuf_append_char(sbuf, '\n');

  return true; // visit all nodes
}


// Get number of links with ltree_get_stats() first
void ltree_write_ctp(LinkTree *tree, const char *kmer, size_t num_links,
                     StrBuf *sbuf)
{
  strbuf_append_str(sbuf, kmer);
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, num_links);
  strbuf_append_char(sbuf, '\n');
  ltree_visit_links(tree, _ltree_write_ctp_link, sbuf);
}

/* Write as dot format */

static inline bool _ltree_dot_node(LinkJunction *l, LinkTree *tree,
                                   uint32_t depth, void *ptr)
{
  (void)depth;
  StrBuf* sbuf = (StrBuf*)ptr;
  const char *seq = ltree_get_seq(tree, l) + (depth == 0 ? tree->kmer_size : 0);
  size_t i, seqlen = strlen(seq);

  strbuf_sprintf(sbuf, "  node%i [label=\"%s\"]\n", l->id, seqlen ? seq : ".");

  for(i = 0; i < 4; i++) {
    if(l->children[i] < 0 && l->counts[i] > 0) { // child but no junction (leaf node)
      strbuf_sprintf(sbuf, "  node%i%c [label=\"%c\"]\n",
                           l->id, "acgt"[i], "ACGT"[i]);
    }
  }

  return true; // print every node
}

static inline bool _ltree_dot_edges(LinkJunction *l, LinkTree *tree,
                                          uint32_t depth, void *ptr)
{
  (void)tree; (void)depth;
  StrBuf* sbuf = (StrBuf*)ptr;
  size_t i;
  for(i = 0; i < 4; i++) {
    if(l->children[i] >= 0 || l->counts[i] > 0) {
      strbuf_sprintf(sbuf, "  node%i -> node%i%c [label=\" %c %zu\"]\n",
                           l->id,
                           l->children[i] < 0 ? l->id : l->children[i],
                           l->children[i] < 0 ? "acgt"[i] : ' ',
                           "ACGT"[i], l->counts[i]);
    }
  }
  return true; // visit every node
}

static void _ltree_dot_root_kmer(LinkTree *tree, StrBuf *sbuf, bool fw)
{
  LTreeID rootid = fw ? tree->fw_id : tree->rv_id;
  if(rootid >= 0) {
    LinkJunction *root = ltree_get_node(tree, rootid);
    strbuf_append_str(sbuf, "  kmer_");
    strbuf_append_str(sbuf, fw ? "fw" : "rv");
    strbuf_append_str(sbuf, "[label=\"");
    strbuf_append_strn(sbuf, ltree_get_seq(tree,root), tree->kmer_size);
    strbuf_append_str(sbuf, fw ? " (F)" : " (R)");
    strbuf_append_str(sbuf, "\"]\n");
  }
}

void ltree_write_dot(LinkTree *tree, StrBuf *sbuf)
{
  strbuf_append_str(sbuf, "digraph G {\n");
  strbuf_append_str(sbuf, "  node [shape=none fontname=\"Courier New\" fontsize=9]\n");
  strbuf_append_str(sbuf, "  edge [shape=none fontname=\"Courier New\" fontsize=9]\n");
  _ltree_dot_root_kmer(tree, sbuf, true);
  _ltree_dot_root_kmer(tree, sbuf, false);
  ltree_visit_nodes(tree, _ltree_dot_node, sbuf);
  if(tree->fw_id >= 0) strbuf_sprintf(sbuf, "  kmer_fw -> node%i\n", tree->fw_id);
  if(tree->rv_id >= 0) strbuf_sprintf(sbuf, "  kmer_rv -> node%i\n", tree->rv_id);
  ltree_visit_nodes(tree, _ltree_dot_edges, sbuf);
  strbuf_append_str(sbuf, "}\n");
}

/* Fetch histogram of coverage */

typedef struct {
  uint64_t *hists; // cast to [col][dist][covg]
  size_t distsize, covgsize;
} LTreeHists;

static inline bool _ltree_update_covg_hists(LinkJunction *l, uint8_t base,
                                            LinkTree *tree,
                                            uint32_t depth, void *ptr)
{
  (void)tree; (void)depth;
  LTreeHists *d = (LTreeHists*)ptr;
  size_t covg = l->counts[base];
  if(l->dist >= d->distsize) return false;
  uint64_t (*hists)[d->covgsize] = (uint64_t (*)[d->covgsize])d->hists;
  hists[l->dist][MIN2(covg,d->covgsize-1)]++;
  return true; // Visit all nodes
}

/**
 * @param hists should be of size distsize * covgsize
 */
void ltree_update_covg_hists(LinkTree *tree, uint64_t *hists,
                             size_t distsize, size_t covgsize)
{
  ctx_assert(hists && distsize > 0 && covgsize > 0);
  LTreeHists data = {.hists = hists, .distsize = distsize, .covgsize = covgsize};
  ltree_visit_links(tree, _ltree_update_covg_hists, &data);
}
