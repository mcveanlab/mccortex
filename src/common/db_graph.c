#include "global.h"
#include "util.h"
#include "binary_kmer.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "path_store.h"
#include "graph_format.h"

dBGraph* db_graph_alloc(dBGraph *db_graph, size_t kmer_size,
                        size_t num_of_cols, size_t num_edge_cols,
                        uint64_t capacity)
{
  size_t i;
  dBGraph tmp = {.kmer_size = kmer_size, .num_of_cols = num_of_cols,
                 .num_edge_cols = num_edge_cols, .num_of_cols_used = 0,
                 .col_edges = NULL, .col_covgs = NULL,
                 .node_in_cols = NULL, .kmer_paths = NULL, .readstrt = NULL};

  memcpy(db_graph, &tmp, sizeof(dBGraph));
  hash_table_alloc(&db_graph->ht, capacity);

  db_graph->ginfo = malloc2(num_of_cols * sizeof(GraphInfo));
  for(i = 0; i < num_of_cols; i++)
    graph_info_alloc(db_graph->ginfo + i);

  char capacity_str[100];
  ulong_to_str(capacity, capacity_str);

  status("[graph] kmer-size: %zu; colours: %zu; capacity: %s\n",
         kmer_size, num_of_cols, capacity_str);

  return db_graph;
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
// Pruning
//

static inline void prune_node_without_edges(dBGraph *db_graph, hkey_t node)
{
  if(db_graph->col_covgs != NULL)
    db_node_zero_covgs(db_graph, node);

  if(db_graph->node_in_cols != NULL) {
    Colour col;
    for(col = 0; col < db_graph->num_of_cols; col++)
      db_node_del_col(db_graph, node, col);
  }

  hash_table_delete(&db_graph->ht, node);
}

static void prune_nodes_lacking_flag(hkey_t node, dBGraph *db_graph,
                                     uint64_t *flags)
{
  assert(db_graph->num_of_cols == 1 && db_graph->col_edges != NULL);

  if(bitset_has(flags, node))
  {
    // Check edges
    Orientation orient, next_orient;
    Nucleotide nuc;
    hkey_t next_node;

    Edges keep_edges = db_graph->col_edges[node];
    BinaryKmer bkmer = db_node_bkmer(db_graph, node);

    for(orient = 0; orient < 2; orient++)
    {
      for(nuc = 0; nuc < 4; nuc++)
      {
        if(edges_has_edge(keep_edges, nuc, orient))
        {
          db_graph_next_node(db_graph, bkmer, nuc, orient,
                             &next_node, &next_orient);

          if(!bitset_has(flags, next_node))
          {
            // Next node fails filter - remove edge
            keep_edges = edges_del_edge(keep_edges, nuc, orient);
          }
        }
      }
    }

    db_graph->col_edges[node] = keep_edges;
  }
  else
  {
    db_graph->col_edges[node] = 0;
    prune_node_without_edges(db_graph, node);
  }
}

static void prune_nodes_lacking_flag_no_edges(hkey_t node, dBGraph *db_graph,
                                              uint64_t *flags)
{
  if(!bitset_has(flags, node))
    prune_node_without_edges(db_graph, node);
}

// If element_has_flag(node, flag) is not true, reset flags
// Remove edges to nodes where !db_node_has_flag(node, flag)
void db_graph_prune_nodes_lacking_flag(dBGraph *db_graph, uint64_t *flags)
{
  if(db_graph->col_edges != NULL) {
    HASH_TRAVERSE(&db_graph->ht, prune_nodes_lacking_flag, db_graph, flags);
  }
  else {
    HASH_TRAVERSE(&db_graph->ht, prune_nodes_lacking_flag_no_edges,
                  db_graph, flags);
  }
}

// Removed edges from nodes that connected to the given `node`
static void prune_connected_nodes(dBGraph *db_graph, hkey_t node, Edges edges)
{
  BinaryKmer bkmer = db_node_bkmer(db_graph, node);
  hkey_t next_node;
  Orientation or, next_or;
  Nucleotide nuc, lost_nuc;
  Edges remove_edge_mask;
  size_t col;

  for(or = 0; or < 2; or++)
  {
    // Remove edge from: next_node:!next_or -> node:!or
    // when working backwards, or = !or, next_or = !next_or
    // so this is actually if(or == REVERSE)

    if(or == FORWARD) {
      lost_nuc = binary_kmer_first_nuc(bkmer, db_graph->kmer_size);
      lost_nuc = binary_nuc_complement(lost_nuc);
    }
    else lost_nuc = binary_kmer_last_nuc(bkmer);

    for(nuc = 0; nuc < 4; nuc++)
    {
      if(edges_has_edge(edges, nuc, or))
      {
        db_graph_next_node(db_graph, bkmer, nuc, or, &next_node, &next_or);
        remove_edge_mask = nuc_orient_to_edge(lost_nuc, rev_orient(next_or));
        assert((db_node_col_edges_union(db_graph, next_node) & remove_edge_mask)
                 == remove_edge_mask);
        for(col = 0; col < db_graph->num_edge_cols; col++)
          db_node_col_edges(db_graph, col, next_node) &= ~remove_edge_mask;
      }
    }
  }
}

void db_graph_prune_node(dBGraph *db_graph, hkey_t node)
{
  size_t col;
  Edges uedges = db_node_col_edges_union(db_graph, node);
  prune_connected_nodes(db_graph, node, uedges);

  for(col = 0; col < db_graph->num_edge_cols; col++)
    db_node_col_edges(db_graph,col,node) = 0;

  prune_node_without_edges(db_graph, node);
}

void db_graph_prune_supernode(dBGraph *db_graph, hkey_t *nodes, size_t len)
{
  if(len == 0) return;
  if(len == 1) db_graph_prune_node(db_graph, nodes[0]);
  else {
    db_graph_prune_node(db_graph, nodes[0]);
    db_graph_prune_node(db_graph, nodes[len-1]);

    size_t i, col;
    for(i = 1; i+1 < len; i++)
      prune_node_without_edges(db_graph, nodes[i]);

    for(i = 1; i+1 < len; i++)
      for(col = 0; col < db_graph->num_edge_cols; col++)
        db_node_col_edges(db_graph, col, nodes[i]) = 0;
  }
}

//
// Functions applying to whole graph
//

// Remove nodes that are in the hash table but not assigned any colours
static void db_graph_remove_node_if_uncoloured(hkey_t node, dBGraph *db_graph)
{
  Colour col = 0;
  while(col < db_graph->num_of_cols && !db_node_has_col(db_graph, node, col))
    col++;

  if(col == db_graph->num_of_cols)
    db_graph_prune_node(db_graph, node);
}

void db_graph_remove_uncoloured_nodes(dBGraph *db_graph)
{
  HASH_TRAVERSE(&db_graph->ht, db_graph_remove_node_if_uncoloured, db_graph);
}

void db_graph_wipe_colour(dBGraph *db_graph, Colour col)
{
  size_t i, capacity = db_graph->ht.capacity, num_of_cols = db_graph->num_of_cols;

  if(db_graph->node_in_cols != NULL)
  {
    size_t words = round_bits_to_words64(capacity);
    for(i = 0; i < words; i++)
      db_graph->node_in_cols[num_of_cols*i+col] = 0;
  }

  Edges (*col_edges)[num_of_cols] = (Edges (*)[num_of_cols])db_graph->col_edges;
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
  Edges edge, *edges = &db_node_col_edges(db_graph,node,0), iedges = edges[0];
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
