#include "global.h"
#include "util.h"
#include "binary_kmer.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "path_store.h"
#include "packed_path.h"
#include "graph_format.h"

static void db_graph_status(const dBGraph *db_graph)
{
  char capacity_str[100];
  ulong_to_str(db_graph->ht.capacity, capacity_str);
  status("[graph] kmer-size: %zu; colours: %zu; capacity: %s\n",
         db_graph->kmer_size, db_graph->num_of_cols, capacity_str);
}

void db_graph_alloc(dBGraph *db_graph, size_t kmer_size,
                    size_t num_of_cols, size_t num_edge_cols,
                    uint64_t capacity)
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
                 .kmer_paths = NULL,
                 .path_kmer_locks = NULL,
                 .readstrt = NULL};

  ctx_assert2(kmer_size >= MIN_KMER_SIZE, "kmer size: %zu", kmer_size);
  ctx_assert2(kmer_size <= MAX_KMER_SIZE, "kmer size: %zu", kmer_size);

  hash_table_alloc(&tmp.ht, capacity);
  memset(&tmp.pdata, 0, sizeof(PathStore));

  tmp.ginfo = malloc2(num_of_cols * sizeof(GraphInfo));
  for(i = 0; i < num_of_cols; i++)
    graph_info_alloc(tmp.ginfo + i);

  memcpy(db_graph, &tmp, sizeof(dBGraph));
  db_graph_status(db_graph);
}

void db_graph_realloc(dBGraph *graph, size_t num_of_cols, size_t num_edge_cols)
{
  // Adjust ginfo size
  size_t i;
  if(graph->num_of_cols == num_of_cols) return;
  if(graph->num_of_cols < num_of_cols) { // Grow
    graph->ginfo = realloc2(graph->ginfo, num_of_cols*sizeof(GraphInfo));
    for(i = graph->num_of_cols; i < num_of_cols; i++)
      graph_info_alloc(graph->ginfo + i);
  }
  else if(graph->num_of_cols > num_of_cols) { // Shrink
    for(i = num_of_cols; i < graph->num_of_cols; i++)
      graph_info_dealloc(graph->ginfo + i);
  }

  dBGraph tmp = {.ht = graph->ht,
                 .kmer_size = graph->kmer_size,
                 .num_of_cols = num_of_cols,
                 .num_edge_cols = num_edge_cols,
                 .num_of_cols_used = graph->num_of_cols_used,
                 .bktlocks = graph->bktlocks,
                 .ginfo = graph->ginfo,
                 .col_edges = graph->col_edges,
                 .col_covgs = graph->col_covgs,
                 .node_in_cols = graph->node_in_cols,
                 .pdata = graph->pdata,
                 .kmer_paths = graph->kmer_paths,
                 .path_kmer_locks = graph->path_kmer_locks,
                 .readstrt = graph->readstrt};

  memcpy(graph, &tmp, sizeof(dBGraph));
  db_graph_status(graph);
}

void db_graph_dealloc(dBGraph *db_graph)
{
  size_t i;
  hash_table_dealloc(&db_graph->ht);
  for(i = 0; i < db_graph->num_of_cols; i++)
    graph_info_dealloc(db_graph->ginfo+i);
  free(db_graph->ginfo);
}

//
// Add to the de bruijn graph
//

void db_graph_update_node_mt(dBGraph *db_graph, dBNode node, Colour col)
{
  if(db_graph->node_in_cols != NULL) db_node_set_col_mt(db_graph, node.key, col);
  if(db_graph->col_covgs != NULL) db_node_increment_coverage_mt(db_graph, node.key, col);
}

// Thread safe
// Note: node may alreay exist in the graph
dBNode db_graph_find_or_add_node_mt(dBGraph *db_graph, BinaryKmer bkmer,
                                    boolean *found)
{
  BinaryKmer bkey;
  hkey_t hkey;
  Orientation orient;

  bkey = bkmer_get_key(bkmer, db_graph->kmer_size);
  hkey = hash_table_find_or_insert_mt(&db_graph->ht, bkey, found,
                                      db_graph->bktlocks);
  orient = bkmer_get_orientation(bkey, bkmer);

  dBNode node = {.key = hkey, .orient = orient};
  return node;
}

dBNode db_graph_find(const dBGraph *db_graph, BinaryKmer bkmer)
{
  dBNode node;
  BinaryKmer bkey;
  bkey = bkmer_get_key(bkmer, db_graph->kmer_size);
  node.key = hash_table_find(&db_graph->ht, bkey);
  node.orient = bkmer_get_orientation(bkmer, bkey);
  return node;
}

// Thread safe
// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge_mt(dBGraph *db_graph, Colour col, dBNode src, dBNode tgt)
{
  if(db_graph->col_edges == NULL) return;

  Nucleotide lhs_nuc, rhs_nuc, lhs_nuc_rev;
  lhs_nuc = db_node_get_first_nuc(src, db_graph);
  rhs_nuc = db_node_get_last_nuc(tgt, db_graph);

  lhs_nuc_rev = dna_nuc_complement(lhs_nuc);

  db_node_set_col_edge_mt(db_graph, src.key, col, rhs_nuc,      src.orient);
  db_node_set_col_edge_mt(db_graph, tgt.key, col, lhs_nuc_rev, !tgt.orient);
}

// For debugging + healthcheck
boolean db_graph_check_edges(const dBGraph *db_graph, dBNode src, dBNode tgt)
{
  Nucleotide lhs_nuc, rhs_nuc, lhs_nuc_rev;
  lhs_nuc = db_node_get_first_nuc(src, db_graph);
  rhs_nuc = db_node_get_last_nuc(tgt, db_graph);

  lhs_nuc_rev = dna_nuc_complement(lhs_nuc);

  Edges src_uedges = db_node_get_edges_union(db_graph, src.key);
  Edges tgt_uedges = db_node_get_edges_union(db_graph, tgt.key);

  return edges_has_edge(src_uedges, rhs_nuc,      src.orient) &&
         edges_has_edge(tgt_uedges, lhs_nuc_rev, !tgt.orient);
}

//
// Graph Traversal
//

dBNode db_graph_next_node(const dBGraph *db_graph,
                          const BinaryKmer node_bkey, Nucleotide next_nuc,
                          Orientation orient)
{
  size_t kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer;

  if(orient == FORWARD)
    bkmer = binary_kmer_left_shift_add(node_bkey, kmer_size, next_nuc);
  else {
    next_nuc = dna_nuc_complement(next_nuc);
    bkmer = binary_kmer_right_shift_add(node_bkey, kmer_size, next_nuc);
  }

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

static inline boolean check_node(hkey_t node, const dBGraph *db_graph)
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

  // Check all edges are reciprical
  for(i = 0; i < nfw_edges; i++)
    if(!db_graph_check_edges(db_graph, fwnode, fwnodes[i])) return false;
  for(i = 0; i < nrv_edges; i++)
    if(!db_graph_check_edges(db_graph, rvnode, rvnodes[i])) return false;

  return true;
}

void db_graph_healthcheck(const dBGraph *db_graph)
{
  HASH_ITERATE(&db_graph->ht, check_node, db_graph);
}

//
// Functions applying to whole graph
//

// BEWARE: if num_edge_cols == 1, edges in all colours will be effectively wiped
void db_graph_wipe_colour(dBGraph *db_graph, Colour col)
{
  status("Wiping graph colour %zu", (size_t)col);

  Edges (*col_edges)[db_graph->num_edge_cols];
  Covg (*col_covgs)[db_graph->num_of_cols];
  const size_t capacity = db_graph->ht.capacity;
  size_t i;

  graph_info_init(&db_graph->ginfo[0]);

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
  size_t col, kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer, bkey, node_bkey = db_node_get_bkmer(db_graph, node);
  Orientation orient;
  Nucleotide nuc;
  hkey_t next;
  Edges edge, *edges = &db_node_edges(db_graph,node,0), iedges = edges[0];
  boolean node_has_col[db_graph->num_edge_cols];

  for(col = 0; col < db_graph->num_edge_cols; col++) {
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

        bkey = bkmer_get_key(bkmer, kmer_size);
        next = hash_table_find(&db_graph->ht, bkey);

        if(next != HASH_NOT_FOUND)
          for(col = 0; col < db_graph->num_edge_cols; col++)
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

//
// Kmer paths
//

void db_graph_dump_paths_by_kmer(const dBGraph *db_graph)
{
  const PathStore *pstore = &db_graph->pdata;
  size_t kmer_size = db_graph->kmer_size;
  char str[MAX_KMER_SIZE+1];
  hkey_t hkey;
  PathIndex pindex;
  PathLen plen;
  Orientation orient, porient;
  boolean first;
  const uint8_t *path;

  printf("\n-------- paths --------\n");

  for(hkey = 0; hkey < db_graph->ht.capacity; hkey++) {
    if(db_graph_node_assigned(db_graph, hkey)) {
      binary_kmer_to_str(db_node_get_bkmer(db_graph, hkey), kmer_size, str);
      for(orient = 0; orient < 2; orient++) {
        pindex = db_node_paths(db_graph, hkey);
        first = true;
        while(pindex != PATH_NULL) {
          path = pstore->store+pindex;
          packedpath_get_len_orient(path, pstore->colset_bytes, &plen, &porient);
          if(porient == orient) {
            if(first) { printf("%s:%i\n", str, orient); first = false; }
            path_store_print_path(pstore, pindex);
          }
          pindex = packedpath_get_prev(pstore->store+pindex);
        }
      }
    }
  }

  printf("-----------------------\n\n");
}

// call seed_random() before any calls to this function please
hkey_t db_graph_rand_node(const dBGraph *db_graph)
{
  uint64_t capacity = db_graph->ht.capacity;
  BinaryKmer *table = db_graph->ht.table;
  hkey_t hkey;

  if(capacity == 0) {
    warn("No entries in hash table - cannot select random");
    return HASH_NOT_FOUND;
  }

  while(1)
  {
    hkey = (hkey_t)((rand() / (double)RAND_MAX) * capacity);
    if(HASH_ENTRY_ASSIGNED(table[hkey])) return hkey;
  }
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
