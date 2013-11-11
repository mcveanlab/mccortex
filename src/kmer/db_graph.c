#include "global.h"
#include "util.h"
#include "binary_kmer.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "path_store.h"
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
  dBGraph tmp = {.kmer_size = kmer_size, .num_of_cols = num_of_cols,
                 .num_edge_cols = num_edge_cols, .num_of_cols_used = 0,
                 .ginfo = NULL, .col_edges = NULL, .col_covgs = NULL,
                 .node_in_cols = NULL, .kmer_paths = NULL, .readstrt = NULL};

  hash_table_alloc(&tmp.ht, capacity);

  tmp.ginfo = malloc2(num_of_cols * sizeof(GraphInfo));
  for(i = 0; i < num_of_cols; i++)
    graph_info_alloc(tmp.ginfo + i);

  memcpy(db_graph, &tmp, sizeof(dBGraph));
  db_graph_status(db_graph);
}

void db_graph_realloc(dBGraph *graph, size_t num_of_cols, size_t num_edge_cols)
{
  // Fix ginfo
  size_t i;
  if(graph->num_of_cols < num_of_cols) { // Grow
    graph->ginfo = realloc2(graph->ginfo, num_of_cols*sizeof(GraphInfo));
    for(i = graph->num_of_cols; i < num_of_cols; i++)
      graph_info_alloc(graph->ginfo + i);
  }
  else if(graph->num_of_cols > num_of_cols) { // Shrink
    for(i = num_of_cols; i < graph->num_of_cols; i++)
      graph_info_dealloc(graph->ginfo + i);
  }
  else return; // Same size

  dBGraph tmp = {.ht = graph->ht, .kmer_size = graph->kmer_size,
                 .num_of_cols = num_of_cols, .num_edge_cols = num_edge_cols,
                 .num_of_cols_used = graph->num_of_cols_used,
                 .ginfo = graph->ginfo, .col_edges = graph->col_edges,
                 .col_covgs = graph->col_covgs, .node_in_cols = graph->node_in_cols,
                 .kmer_paths = graph->kmer_paths, .readstrt = graph->readstrt};

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

// Note: node may alreay exist in the graph
hkey_t db_graph_find_or_add_node(dBGraph *db_graph, BinaryKmer bkey, Colour col)
{
  boolean found;
  hkey_t node = hash_table_find_or_insert(&db_graph->ht, bkey, &found);
  if(db_graph->node_in_cols != NULL) db_node_set_col(db_graph, node, col);
  if(db_graph->col_covgs != NULL) db_node_increment_coverage(db_graph, node, col);
  return node;
}

// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge(dBGraph *db_graph, Colour colour,
                       hkey_t src_node, hkey_t tgt_node,
                       Orientation src_orient, Orientation tgt_orient)
{
  if(db_graph->col_edges == NULL) return;

  BinaryKmer src_bkmer = db_node_bkmer(db_graph, src_node);
  BinaryKmer tgt_bkmer = db_node_bkmer(db_graph, tgt_node);

  Nucleotide lhs_nuc, rhs_nuc;
  lhs_nuc = db_node_first_nuc(src_bkmer, src_orient, db_graph->kmer_size);
  rhs_nuc = db_node_last_nuc(tgt_bkmer, tgt_orient, db_graph->kmer_size);

  Nucleotide lhs_nuc_rev = binary_nuc_complement(lhs_nuc);
  Orientation tgt_orient_opp = opposite_orientation(tgt_orient);

  db_node_set_col_edge(db_graph, colour, src_node, rhs_nuc, src_orient);
  db_node_set_col_edge(db_graph, colour, tgt_node, lhs_nuc_rev, tgt_orient_opp);
}

//
// Graph Traversal
//

void db_graph_next_node(const dBGraph *db_graph,
                        BinaryKmer bkmer, Nucleotide next_nuc,
                        Orientation orient,
                        hkey_t *next_node, Orientation *next_orient)
{
  size_t kmer_size = db_graph->kmer_size;
  if(orient == FORWARD) binary_kmer_left_shift_add(&bkmer, kmer_size, next_nuc);
  else binary_kmer_right_shift_add(&bkmer, kmer_size, binary_nuc_complement(next_nuc));
  BinaryKmer bkey = db_node_get_key(bkmer, kmer_size);
  *next_node = hash_table_find(&db_graph->ht, bkey);
  *next_orient = db_node_get_orientation(bkmer, bkey) ^ orient;
  assert(*next_node != HASH_NOT_FOUND);
}

size_t db_graph_next_nodes(const dBGraph *db_graph,
                           BinaryKmer bkmer, Orientation orient, Edges edges,
                           hkey_t nodes[4], Orientation orients[4],
                           Nucleotide fw_nucs[4])
{
  size_t count = 0, kmer_size = db_graph->kmer_size;
  Edges tmp_edge;
  Nucleotide nuc;
  BinaryKmer bkey;

  edges = edges_with_orientation(edges, orient);
  if(orient == FORWARD) binary_kmer_left_shift_one_base(&bkmer, kmer_size);
  else binary_kmer_right_shift_one_base(&bkmer);

  for(tmp_edge = 0x1, nuc = 0; nuc < 4; tmp_edge <<= 1, nuc++) {
    if(edges & tmp_edge) {
      if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
      else binary_kmer_set_first_nuc(&bkmer, binary_nuc_complement(nuc), kmer_size);
      bkey = db_node_get_key(bkmer, kmer_size);
      nodes[count] = hash_table_find(&db_graph->ht, bkey);
      orients[count] = db_node_get_orientation(bkmer, bkey) ^ orient;
      fw_nucs[count] = nuc;
      count++;
    }
  }

  return count;
}


//
// Functions applying to whole graph
//

void db_graph_wipe_colour(dBGraph *db_graph, Colour col)
{
  size_t i, capacity, num_of_cols, num_edge_cols;

  capacity = db_graph->ht.capacity;
  num_of_cols = db_graph->num_of_cols;
  num_edge_cols = db_graph->num_edge_cols;

  if(db_graph->node_in_cols != NULL)
  {
    size_t words = round_bits_to_words64(capacity);
    for(i = 0; i < words; i++)
      db_graph->node_in_cols[num_of_cols*i+col] = 0;
  }

  Edges (*col_edges)[num_edge_cols] = (Edges (*)[num_edge_cols])db_graph->col_edges;
  Covg (*col_covgs)[num_of_cols] = (Covg (*)[num_of_cols])db_graph->col_covgs;

  if(db_graph->col_covgs != NULL) {
    for(i = 0; i < capacity; i++)
      col_covgs[i][col] = 0;
  }

  if(db_graph->col_edges != NULL) {
    for(i = 0; i < capacity; i++)
      col_edges[i][col] = 0;
  }
}

static inline void add_all_edges(hkey_t node, dBGraph *db_graph)
{
  size_t col, kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer, bkey, node_bkey = db_node_bkmer(db_graph, node);
  Orientation orient;
  Nucleotide nuc;
  hkey_t next;
  Edges edge, *edges = &db_node_col_edges(db_graph,0,node), iedges = edges[0];
  boolean node_has_col[db_graph->num_edge_cols];

  for(col = 0; col < db_graph->num_edge_cols; col++) {
    iedges &= edges[col];
    node_has_col[col] = db_node_has_col(db_graph, node, col);
  }

  for(orient = 0; orient < 2; orient++)
  {
    bkmer = node_bkey;
    if(orient == FORWARD) binary_kmer_left_shift_one_base(&bkmer, kmer_size);
    else binary_kmer_right_shift_one_base(&bkmer);

    for(nuc = 0; nuc < 4; nuc++)
    {
      edge = nuc_orient_to_edge(nuc, orient);

      // Check edge is not is all colours
      if(!(edge & iedges))
      {
        if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
        else binary_kmer_set_first_nuc(&bkmer, binary_nuc_complement(nuc), kmer_size);

        bkey = db_node_get_key(bkmer, kmer_size);
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
  assert(db_graph->num_of_cols == db_graph->num_edge_cols);
  HASH_TRAVERSE(&db_graph->ht, add_all_edges, db_graph);
}

//
// Kmer paths
//

void db_graph_dump_paths_by_kmer(const dBGraph *db_graph)
{
  const PathStore *paths = &db_graph->pdata;
  uint32_t kmer_size = db_graph->kmer_size;
  char str[MAX_KMER_SIZE+1];
  hkey_t node;
  PathIndex index, prev_index;
  PathLen len;
  Orientation orient, porient;
  boolean first;

  printf("\n-------- paths --------\n");

  for(node = 0; node < db_graph->ht.capacity; node++) {
    if(db_graph_node_assigned(db_graph, node)) {
      binary_kmer_to_str(db_node_bkmer(db_graph, node), kmer_size, str);
      for(orient = 0; orient < 2; orient++) {
        index = db_node_paths(db_graph, node);
        first = true;
        while(index != PATH_NULL) {
          prev_index = path_store_prev(paths, index);
          path_store_len_orient(paths, index, &len, &porient);
          if(porient == orient) {
            if(first) { printf("%s:%i\n", str, orient); first = false; }
            path_store_print_path(paths, index);
          }
          index = prev_index;
        }
      }
    }
  }

  printf("-----------------------\n\n");
}

hkey_t db_graph_rand_node(const dBGraph *db_graph)
{
  uint64_t capacity = db_graph->ht.capacity;
  BinaryKmer *table = db_graph->ht.table;
  hkey_t node;

  if(capacity == 0) {
    warn("No entries in hash table - cannot select random");
    return HASH_NOT_FOUND;
  }

  while(1)
  {
    node = ((double)rand() / RAND_MAX) * capacity;
    if(HASH_ENTRY_ASSIGNED(table[node])) return node;
  }
  // while(node < capacity && !HASH_ENTRY_ASSIGNED(table[node])) node++;
  // if(node == capacity) {
  //   node = 0;
  //   while(node < capacity && !HASH_ENTRY_ASSIGNED(table[node])) node++;
  // }
  return node;
}
