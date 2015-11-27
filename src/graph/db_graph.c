#include "global.h"
#include "util.h"
#include "binary_kmer.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"

static void db_graph_status(const dBGraph *db_graph)
{
  char capacity_str[100];
  ulong_to_str(db_graph->ht.capacity, capacity_str);
  status("[graph] kmer-size: %zu; colours: %zu; capacity: %s\n",
         db_graph->kmer_size, db_graph->num_of_cols, capacity_str);
}

const int DBG_ALLOC_EDGES       =  1;
const int DBG_ALLOC_COVGS       =  2;
const int DBG_ALLOC_BKTLOCKS    =  4;
const int DBG_ALLOC_READSTRT    =  8;
const int DBG_ALLOC_NODE_IN_COL = 16;

// alloc_flags specifies where fields to malloc. OR together DBG_ALLOC_* values
void db_graph_alloc(dBGraph *db_graph, size_t kmer_size,
                    size_t num_of_cols, size_t num_edge_cols,
                    uint64_t capacity, int alloc_flags)
{
  size_t i;
  dBGraph tmp = {.kmer_size = kmer_size,
                 .num_of_cols = num_of_cols,
                 .num_edge_cols = num_edge_cols,
                 .num_of_cols_used = 0,
                 .bktlocks = NULL,
                 .ginfo = NULL,
                 .col_edges = NULL,
                 .col_covgs = NULL,
                 .node_in_cols = NULL,
                 .readstrt = NULL};

  ctx_assert(num_of_cols > 0);
  ctx_assert(num_edge_cols == 0 || num_edge_cols == 1 || num_edge_cols == num_of_cols);
  ctx_assert(capacity > 0);
  ctx_assert2(kmer_size >= MIN_KMER_SIZE, "kmer size: %zu", kmer_size);
  ctx_assert2(kmer_size <= MAX_KMER_SIZE, "kmer size: %zu", kmer_size);

  hash_table_alloc(&tmp.ht, capacity);
  memset(&tmp.gpstore, 0, sizeof(GPathStore));

  tmp.ginfo = ctx_calloc(num_of_cols, sizeof(GraphInfo));
  for(i = 0; i < num_of_cols; i++)
    graph_info_alloc(&tmp.ginfo[i]);

  if(alloc_flags & DBG_ALLOC_EDGES)
    tmp.col_edges = ctx_calloc(tmp.ht.capacity * num_edge_cols, sizeof(Edges));

  if(alloc_flags & DBG_ALLOC_COVGS)
    tmp.col_covgs = ctx_calloc(tmp.ht.capacity * num_of_cols, sizeof(Covg));

  if(alloc_flags & DBG_ALLOC_BKTLOCKS)
    tmp.bktlocks = ctx_calloc(roundup_bits2bytes(tmp.ht.num_of_buckets), 1);

  // 1 bit for forward, 1 bit for reverse per kmer
  if(alloc_flags & DBG_ALLOC_READSTRT)
    tmp.readstrt = ctx_calloc(roundup_bits2bytes(tmp.ht.capacity)*2, 1);

  if(alloc_flags & DBG_ALLOC_NODE_IN_COL) {
    size_t bytes_per_col = roundup_bits2bytes(tmp.ht.capacity);
    tmp.node_in_cols = ctx_calloc(bytes_per_col*num_of_cols, 1);
  }

  memcpy(db_graph, &tmp, sizeof(dBGraph));
  db_graph_status(db_graph);
}

// Free memory used by all fields as well
void db_graph_dealloc(dBGraph *db_graph)
{
  size_t i;

  hash_table_dealloc(&db_graph->ht);

  for(i = 0; i < db_graph->num_of_cols; i++)
    graph_info_dealloc(db_graph->ginfo+i);
  ctx_free(db_graph->ginfo);

  ctx_free(db_graph->bktlocks);
  ctx_free(db_graph->col_covgs); // num_of_cols * capacity
  ctx_free(db_graph->col_edges); // num_col_edges * capacity
  ctx_free(db_graph->node_in_cols);
  ctx_free(db_graph->readstrt);

  gpath_hash_dealloc(&db_graph->gphash);
  gpath_store_dealloc(&db_graph->gpstore);

  memset(db_graph, 0, sizeof(dBGraph));
}

//
// Add to the de bruijn graph
//

void db_graph_update_node_mt(dBGraph *db_graph, dBNode node, Colour col)
{
  if(db_graph->node_in_cols != NULL) db_node_set_col_mt(db_graph, node.key, col);
  if(db_graph->col_covgs != NULL) db_node_increment_coverage_mt(db_graph, node.key, col);
}

// Not thread safe, use db_graph_find_or_add_node_mt for that
// Note: node may alreay exist in the graph
dBNode db_graph_find_or_add_node(dBGraph *db_graph, BinaryKmer bkmer,
                                 bool *foundptr)
{
  BinaryKmer bkey = binary_kmer_get_key(bkmer, db_graph->kmer_size);
  hkey_t hkey = hash_table_find_or_insert(&db_graph->ht, bkey, foundptr);
  return (dBNode){.key = hkey, .orient = bkmer_get_orientation(bkey, bkmer)};
}

// Thread safe
// Note: node may alreay exist in the graph
dBNode db_graph_find_node_mt(dBGraph *db_graph, BinaryKmer bkmer)
{
  BinaryKmer bkey = binary_kmer_get_key(bkmer, db_graph->kmer_size);
  hkey_t hkey = hash_table_find_mt(&db_graph->ht, bkey, db_graph->bktlocks);
  return (dBNode){.key = hkey, .orient = bkmer_get_orientation(bkey, bkmer)};
}

dBNode db_graph_find_or_add_node_mt(dBGraph *db_graph, BinaryKmer bkmer,
                                    bool *foundptr)
{
  BinaryKmer bkey = binary_kmer_get_key(bkmer, db_graph->kmer_size);
  hkey_t hkey = hash_table_find_or_insert_mt(&db_graph->ht, bkey, foundptr,
                                             db_graph->bktlocks);

  return (dBNode){.key = hkey, .orient = bkmer_get_orientation(bkey, bkmer)};
}

dBNode db_graph_find_str(const dBGraph *db_graph, const char *str)
{
  BinaryKmer bkmer;
  bkmer = binary_kmer_from_str(str, db_graph->kmer_size);
  return db_graph_find(db_graph, bkmer);
}

dBNode db_graph_find_node(const dBGraph *db_graph, BinaryKmer bkmer)
{
  BinaryKmer bkey = binary_kmer_get_key(bkmer, db_graph->kmer_size);
  hkey_t hkey = hash_table_find(&db_graph->ht, bkey);
  return (dBNode){.key = hkey, .orient = bkmer_get_orientation(bkey, bkmer)};
}

// Thread safe
// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge_mt(dBGraph *db_graph, Colour col, dBNode src, dBNode tgt)
{
  if(db_graph->col_edges == NULL) return;

  ctx_assert(col < db_graph->num_edge_cols);

  Nucleotide lhs_nuc, rhs_nuc, lhs_nuc_rev;
  lhs_nuc = db_node_get_first_nuc(src, db_graph);
  rhs_nuc = db_node_get_last_nuc(tgt, db_graph);

  lhs_nuc_rev = dna_nuc_complement(lhs_nuc);

  db_node_set_col_edge_mt(db_graph, src.key, col, rhs_nuc,      src.orient);
  db_node_set_col_edge_mt(db_graph, tgt.key, col, lhs_nuc_rev, !tgt.orient);
}

// For debugging + healthcheck
// dies with message if edge from src->tgt and not tgt->src or vice-versa
// returns true if edges would be added by infer edges
bool db_graph_check_edges(const dBGraph *db_graph, dBNode src, dBNode tgt)
{
  Nucleotide lhs_nuc, rhs_nuc, lhs_nuc_rev;
  Edges src_edges, tgt_edges;
  size_t col;
  bool fw_edge, rv_edge, missing_edges = false;

  lhs_nuc = db_node_get_first_nuc(src, db_graph);
  rhs_nuc = db_node_get_last_nuc(tgt, db_graph);

  lhs_nuc_rev = dna_nuc_complement(lhs_nuc);

  for(col = 0; col < db_graph->num_edge_cols; col++) {
    src_edges = db_node_get_edges(db_graph, src.key, col);
    tgt_edges = db_node_get_edges(db_graph, tgt.key, col);

    fw_edge = edges_has_edge(src_edges, rhs_nuc,      src.orient);
    rv_edge = edges_has_edge(tgt_edges, lhs_nuc_rev, !tgt.orient);

    if(fw_edge != rv_edge) die("Missing edge pair");

    missing_edges |= !fw_edge &&
                     db_node_has_col(db_graph, src.key, col) &&
                     db_node_has_col(db_graph, tgt.key, col);
  }

  return missing_edges;
}

// Returns true on missing edges, false otherwise
bool db_graph_check_all_edges(const dBGraph *db_graph,
                              const dBNode *nodes, size_t num_nodes)
{
  size_t i;
  for(i = 0; i+1 < num_nodes; i++)
    if(db_graph_check_edges(db_graph, nodes[i], nodes[i+1]))
      return true;

  return false;
}

//
// Graph Traversal
//

dBNode db_graph_next_node(const dBGraph *db_graph,
                          const BinaryKmer node_bkey, Nucleotide next_nuc,
                          Orientation orient)
{
  size_t kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer = bkmer_shift_add_last_nuc(node_bkey, orient,
                                              kmer_size, next_nuc);

  dBNode next_node = db_graph_find(db_graph, bkmer);
  next_node.orient ^= orient;

  ctx_assert(next_node.key != HASH_NOT_FOUND);
  return next_node;
}

uint8_t db_graph_next_nodes(const dBGraph *db_graph, const BinaryKmer node_bkey,
                            Orientation orient, Edges edges,
                            dBNode nodes[4], Nucleotide fw_nucs[4])
{
  const size_t kmer_size = db_graph->kmer_size;
  Edges tmp_edge;
  Nucleotide nuc;
  BinaryKmer bkmer;
  uint8_t count = 0;

  edges = edges_with_orientation(edges, orient);
  bkmer = (orient == FORWARD ? binary_kmer_left_shift_one_base(node_bkey, kmer_size)
                             : binary_kmer_right_shift_one_base(node_bkey));

  for(tmp_edge = 0x1, nuc = 0; nuc < 4; tmp_edge <<= 1, nuc++) {
    if(edges & tmp_edge) {
      if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
      else binary_kmer_set_first_nuc(&bkmer, dna_nuc_complement(nuc), kmer_size);
      nodes[count] = db_graph_find(db_graph, bkmer);
      nodes[count].orient ^= orient;
      fw_nucs[count] = nuc;
      ctx_assert(nodes[count].key != HASH_NOT_FOUND);
      count++;
    }
  }

  return count;
}

uint8_t db_graph_next_nodes_union(const dBGraph *db_graph, dBNode node,
                                  dBNode nodes[4], Nucleotide fw_nucs[4])
{
  BinaryKmer bkey = db_node_get_bkmer(db_graph, node.key);
  Edges edges = db_node_get_edges_union(db_graph, node.key);
  return db_graph_next_nodes(db_graph, bkey, node.orient, edges, nodes, fw_nucs);
}

/**
 * @param colour if > -1: filter next nodes for those in colour, otherwise all next nodes
 * @param fw_nucs is the nuc you would add when walking forward
 * @return Number of nodes added
 */
uint8_t db_graph_next_nodes_in_col(const dBGraph *db_graph,
                                   dBNode node, int colour,
                                   dBNode nodes[4], Nucleotide fw_nucs[4])
{
  ctx_assert2(colour < 0 ||
              (db_graph->num_of_cols == 1 && colour == 0) ||
              db_graph->num_of_cols == db_graph->num_edge_cols ||
              (db_graph->num_of_cols > 1 && db_graph->num_edge_cols == 1 &&
                (db_graph->node_in_cols || db_graph->col_covgs)),
              "col: %i; cols: %zu edges: %zu node_in_cols: %i col_covgs: %i",
              colour, db_graph->num_of_cols, db_graph->num_edge_cols,
              !!db_graph->node_in_cols, !!db_graph->col_covgs);

  size_t i, j;
  Edges edges;
  BinaryKmer bkey;
  uint8_t count;

  if(colour < 0)
    edges = db_node_get_edges_union(db_graph, node.key);
  else if(db_graph->num_edge_cols < db_graph->num_of_cols)
    edges = db_node_edges(db_graph, node.key, 0);
  else
    edges = db_node_edges(db_graph, node.key, colour);

  bkey = db_node_get_bkmer(db_graph, node.key);

  count = db_graph_next_nodes(db_graph, bkey, node.orient, edges,
                              nodes, fw_nucs);

  // Filter next nodes if needed
  // If we allow kmers to exist and not be in any colour when
  //   num_of_cols=1, num_edge_cols=1
  // then we should comment out this if condition
  if(colour >= 0 && db_graph->num_edge_cols < db_graph->num_of_cols)
  {
    for(i = j = 0; i < count; i++) {
      if(( db_graph->node_in_cols && db_node_has_col(db_graph, nodes[i].key, colour)) ||
         (!db_graph->node_in_cols && db_node_covg(db_graph, nodes[i].key, colour) > 0))
      {
        nodes[j] = nodes[i];
        fw_nucs[j] = fw_nucs[i];
        j++;
      }
    }
    count = j;
  }

  return count;
}

/**
 * Get previous nodes in this colour, ignoring the the node we just came from
 * @param node      Current node
 * @param lost_nuc  The first nucleotide of the previous BinaryKmer - the one
 *                  that was lost when traversing to the current node
 * @param colour    Filter down to nodes only in this colour (if >= 0)
 */
uint8_t db_graph_prev_nodes_with_mask(const dBGraph *db_graph, dBNode node,
                                      Nucleotide lost_nuc, int colour,
                                      dBNode prev_nodes[4],
                                      Nucleotide prev_bases[4])
{
  Edges edges, prev_edge;

  if(colour >= 0 && db_graph->num_edge_cols > 1)
    edges = db_node_get_edges(db_graph, node.key, colour);
  else
    edges = db_node_get_edges_union(db_graph, node.key);

  // Remove edge to kmer we came from
  // Can slim down the number of nodes to look up if we can rule out
  // the node we just came from
  lost_nuc = dna_nuc_complement(lost_nuc);
  prev_edge = nuc_orient_to_edge(lost_nuc, rev_orient(node.orient));

  // Some sanity checks
  ctx_assert(edges & prev_edge);
  edges &= ~prev_edge;

  BinaryKmer bkey = db_node_get_bkmer(db_graph, node.key);

  uint8_t i, j, num_prev;

  num_prev = db_graph_next_nodes(db_graph, bkey,
                                 rev_orient(node.orient), edges,
                                 prev_nodes, prev_bases);

  // If we have the ability, slim down nodes by those in this colour
  if(colour >= 0 && db_graph->node_in_cols != NULL) {
    for(i = j = 0; i < num_prev; i++) {
      if(db_node_has_col(db_graph, prev_nodes[i].key, colour)) {
        prev_nodes[j] = prev_nodes[i];
        prev_bases[j] = prev_bases[i];
        j++;
      }
    }
    num_prev = j;
  }

  // Reverse orientation
  for(i = 0; i < num_prev; i++)
    prev_nodes[i].orient = !prev_nodes[i].orient;

  return num_prev;
}

//
//
//

// Check kmer size of a file
// Used when loading graph and path files
void db_graph_check_kmer_size(size_t kmer_size, const char *path)
{
  const size_t mink = MIN_KMER_SIZE, maxk = MAX_KMER_SIZE;
  if(kmer_size < mink || kmer_size > maxk)
    die("Cannot handle kmer size %zu [%zu-%zu; %s]", kmer_size, mink, maxk, path);
  else if(!(kmer_size & 1))
    die("Kmer size appears to be even! %zu [%s]", kmer_size, path);
}

//
// Health check
//

/*
  @param missing_edges edges would be added by infer edges operation
*/
static inline void check_node(hkey_t node, const dBGraph *db_graph,
                              bool *missing_edges_ptr)
{
  Edges edges = db_node_get_edges_union(db_graph, node);
  BinaryKmer bkmer = db_node_get_bkmer(db_graph, node);
  size_t nfw_edges, nrv_edges, i, j;
  dBNode fwnodes[8], rvnodes[8];
  Nucleotide fwnucs[8], rvnucs[8];

  nfw_edges = db_graph_next_nodes(db_graph, bkmer, FORWARD, edges,
                                  fwnodes, fwnucs);

  nrv_edges = db_graph_next_nodes(db_graph, bkmer, REVERSE, edges,
                                  rvnodes, rvnucs);

  for(i = 0; i < nfw_edges && fwnodes[i].key != HASH_NOT_FOUND; i++);
  for(j = 0; j < nrv_edges && rvnodes[j].key != HASH_NOT_FOUND; j++);

  size_t total_edges = nfw_edges + nrv_edges;

  if((unsigned)__builtin_popcount(edges) != total_edges || i+j != total_edges) {
    char seq[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph->kmer_size, seq);
    die("Excess edges on node: %s [%zu,%zu]", seq, nfw_edges, nrv_edges);
  }

  dBNode fwnode = {.key = node, .orient = FORWARD};
  dBNode rvnode = {.key = node, .orient = REVERSE};

  bool missing_edges = false;

  // Check all edges are reciprical
  for(i = 0; i < nfw_edges; i++)
    missing_edges |= db_graph_check_edges(db_graph, fwnode, fwnodes[i]);

  for(i = 0; i < nrv_edges; i++)
    missing_edges |= db_graph_check_edges(db_graph, rvnode, rvnodes[i]);

  *missing_edges_ptr |= missing_edges;
}

void db_graph_healthcheck(const dBGraph *db_graph)
{
  status("Running graph edge check...");
  ctx_assert(db_graph->col_edges != NULL);
  bool missing_edges = false;
  HASH_ITERATE(&db_graph->ht, check_node, db_graph, &missing_edges);
  if(missing_edges) status("  edges would be added with infer edges");
  if(missing_edges) status("  all edges present");
}

//
// Functions applying to whole graph
//

void db_graph_reset(dBGraph *db_graph)
{
  size_t col, capacity = db_graph->ht.capacity;
  size_t ncols = db_graph->num_of_cols, nedgecols = db_graph->num_edge_cols;

  for(col = 0; col < ncols; col++)
    graph_info_init(&db_graph->ginfo[col]);

  hash_table_empty(&db_graph->ht);
  db_graph->num_of_cols_used = 0;

  if(db_graph->col_edges != NULL)
    memset(db_graph->col_edges, 0, nedgecols * sizeof(Edges) * capacity);
  if(db_graph->col_covgs != NULL)
    memset(db_graph->col_covgs, 0, ncols * sizeof(Covg) * capacity);
  if(db_graph->node_in_cols != NULL)
    memset(db_graph->node_in_cols, 0, roundup_bits2bytes(capacity) * ncols);
  if(db_graph->readstrt != NULL)
    memset(db_graph->readstrt, 0, 2 * roundup_bits2bytes(capacity) * ncols);

  gpath_store_reset(&db_graph->gpstore);
}

//
// Stats: Get kmer coverage in each colour
//

typedef struct {
  const dBGraph *db_graph;
  const size_t nthreads;
  uint64_t *nkmers, *sumcov;
} GetKmerCovg;

bool get_kmer_covg(hkey_t hkey, const dBGraph *db_graph,
                   uint64_t *nkmers, uint64_t *sumcov)
{
  size_t col, covg;
  for(col = 0; col < db_graph->num_of_cols; col++) {
    covg = db_node_get_covg(db_graph, hkey, col);
    if(covg) {
      nkmers[col]++;
      sumcov[col] += covg;
    }
  }
  return false; // keep iterating
}

void get_kmer_covg_thread(void *arg, size_t threadid)
{
  GetKmerCovg *d = (GetKmerCovg*)arg;

  size_t col, ncols = d->db_graph->num_of_cols;
  uint64_t *nkmers = ctx_calloc(ncols, sizeof(uint64_t));
  uint64_t *sumcov = ctx_calloc(ncols, sizeof(uint64_t));

  HASH_ITERATE_PART(&d->db_graph->ht, threadid, d->nthreads,
                    get_kmer_covg, d->db_graph, nkmers, sumcov);

  // Add results to array shared with other threads
  for(col = 0; col < ncols; col++) {
    __sync_fetch_and_add((volatile uint64_t*)&d->nkmers[col], nkmers[col]);
    __sync_fetch_and_add((volatile uint64_t*)&d->sumcov[col], sumcov[col]);
  }

  ctx_free(nkmers);
  ctx_free(sumcov);
}

void db_graph_get_kmer_covg(const dBGraph *db_graph, size_t nthreads,
                            uint64_t *nkmers, uint64_t *sumcov)
{
  GetKmerCovg getcov = {.db_graph = db_graph,
                        .nthreads = nthreads,
                        .nkmers = nkmers,
                        .sumcov = sumcov};

  util_multi_thread(&getcov, nthreads, get_kmer_covg_thread);
}

//
// Wipe
//

// BEWARE: if num_edge_cols == 1, edges in all colours will be effectively wiped
void db_graph_wipe_colour(dBGraph *db_graph, Colour col)
{
  status("Wiping graph colour %zu", (size_t)col);

  Edges (*col_edges)[db_graph->num_edge_cols];
  Covg (*col_covgs)[db_graph->num_of_cols];
  const size_t capacity = db_graph->ht.capacity;
  size_t i;

  graph_info_init(&db_graph->ginfo[col]);

  if(db_graph->node_in_cols != NULL)
  {
    size_t nbytes = roundup_bits2bytes(capacity);
    for(i = 0; i < nbytes; i++)
      db_graph->node_in_cols[db_graph->num_of_cols*i+col] = 0;
  }

  col_edges = (Edges (*)[db_graph->num_edge_cols])db_graph->col_edges;
  col_covgs = (Covg (*)[db_graph->num_of_cols])db_graph->col_covgs;

  if(db_graph->col_covgs != NULL) {
    if(db_graph->num_of_cols == 1) {
      memset(db_graph->col_covgs, 0, capacity * sizeof(Covg));
    } else {
      for(i = 0; i < capacity; i++)
        col_covgs[i][col] = 0;
    }
  }

  if(db_graph->col_edges != NULL) {
    if(db_graph->num_edge_cols == 1) {
      memset(db_graph->col_edges, 0, capacity * sizeof(Edges));
    } else {
      for(i = 0; i < capacity; i++)
        col_edges[i][col] = 0;
    }
  }
}

static inline void add_all_edges(hkey_t node, dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size, edgencols = db_graph->num_edge_cols;
  size_t col;
  BinaryKmer bkmer, bkey, node_bkey = db_node_get_bkmer(db_graph, node);
  Orientation orient;
  Nucleotide nuc;
  hkey_t next;
  Edges edge, *edges = &db_node_edges(db_graph,node,0), iedges = edges[0];
  bool node_has_col[edgencols];

  for(col = 0; col < edgencols; col++) {
    iedges &= edges[col];
    node_has_col[col] = db_node_has_col(db_graph, node, col);
  }

  for(orient = 0; orient < 2; orient++)
  {
    bkmer = (orient == FORWARD ? binary_kmer_left_shift_one_base(node_bkey, kmer_size)
                               : binary_kmer_right_shift_one_base(node_bkey));

    for(nuc = 0; nuc < 4; nuc++)
    {
      edge = nuc_orient_to_edge(nuc, orient);

      // Check edge is not is all colours
      if(!(edge & iedges))
      {
        if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
        else binary_kmer_set_first_nuc(&bkmer, dna_nuc_complement(nuc), kmer_size);

        bkey = binary_kmer_get_key(bkmer, kmer_size);
        next = hash_table_find(&db_graph->ht, bkey);

        if(next != HASH_NOT_FOUND)
          for(col = 0; col < edgencols; col++)
            if(node_has_col[col] && db_node_has_col(db_graph, next, col))
              edges[col] |= edge;
      }
    }
  }
}

void db_graph_add_all_edges(dBGraph *db_graph)
{
  ctx_assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  HASH_ITERATE(&db_graph->ht, add_all_edges, db_graph);
}

static bool wipe_kmer_if_no_covg(hkey_t hkey, size_t threadid, void *arg)
{
  (void)threadid;
  dBGraph *db_graph = (dBGraph*)arg;
  size_t col; Covg covg = 0;
  for(col = 0; col < db_graph->num_of_cols; col++)
    covg |= db_node_get_covg(db_graph, hkey, col);
  if(!covg) hash_table_delete(&db_graph->ht, hkey);
  return false; // keep iterating
}

// remove kmers from the graph if they have no coverage
void db_graph_remove_no_covg_kmers(dBGraph *db_graph, size_t nthreads)
{
  hash_table_iterate(&db_graph->ht, nthreads, wipe_kmer_if_no_covg, db_graph);
}

typedef struct {
  Edges *isec_edges;
  dBGraph *db_graph;
  size_t nthreads;
} IntersectEdgesJob;

static void intersect_edges(void *arg, size_t threadid)
{
  IntersectEdgesJob job = *(IntersectEdgesJob*)arg;
  Edges *edges = job.db_graph->col_edges;
  size_t step, start, end, i, j, col, ncols;
  step = job.db_graph->ht.capacity / job.nthreads;
  start = step * threadid;
  end = threadid+1 == job.nthreads ? job.db_graph->ht.capacity : start + step;
  ncols = job.db_graph->num_of_cols;
  for(i = start, j = i*ncols; i < end; i++)
    for(col = 0; col < ncols; col++, j++)
      edges[j] &= job.isec_edges[i];
}

void db_graph_intersect_edges(dBGraph *db_graph, size_t nthreads, Edges *edges)
{
  IntersectEdgesJob job = {.isec_edges = edges,
                           .db_graph = db_graph,
                           .nthreads = nthreads};
  util_multi_thread(&job, nthreads, intersect_edges);
}

// Get a random node from the graph
// call seed_random() before any calls to this function please
// if ntries > 0 and we fail to find a node will return HASH_NOT_FOUND
hkey_t db_graph_rand_node(const dBGraph *db_graph, size_t ntries)
{
  uint64_t capacity = db_graph->ht.capacity;
  BinaryKmer *table = db_graph->ht.table;
  hkey_t hkey;
  size_t i;

  if(capacity == 0) {
    warn("No entries in hash table - cannot select random");
    return HASH_NOT_FOUND;
  }

  for(i = 0; i < ntries; i++)
  {
    hkey = (hkey_t)((rand() / (double)RAND_MAX) * capacity);
    if(HASH_ENTRY_ASSIGNED(table[hkey])) return hkey;
  }

  return HASH_NOT_FOUND;
}

void db_graph_print_kmer2(BinaryKmer bkmer, Covg *covgs, Edges *edges,
                          size_t num_of_cols, size_t kmer_size, FILE *fout)
{
  char bkmerstr[MAX_KMER_SIZE+1], edgesstr[10];
  size_t i;

  binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
  fputs(bkmerstr, fout);

  // Print covgs
  for(i = 0; i < num_of_cols; i++)
    fprintf(fout, " %u", covgs[i]);

  // Print edges
  for(i = 0; i < num_of_cols; i++) {
    fputc(' ', fout);
    fputs(db_node_get_edges_str(edges[i], edgesstr), fout);
  }

  fputc('\n', fout);
}

void db_graph_print_kmer(hkey_t node, dBGraph *db_graph, FILE *fout)
{
  BinaryKmer bkmer = db_node_get_bkmer(db_graph, node);
  Covg *covgs = &db_node_covg(db_graph, node, 0);
  Edges *edges = &db_node_edges(db_graph, node, 0);

  db_graph_print_kmer2(bkmer, covgs, edges,
                       db_graph->num_of_cols, db_graph->kmer_size,
                       fout);
}
