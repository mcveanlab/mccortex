#include "global.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "binary_paths.h"

dBGraph* db_graph_alloc(dBGraph *db_graph, uint32_t kmer_size, uint64_t capacity)
{
  memset(db_graph, 0, sizeof(dBGraph));
  db_graph->kmer_size = kmer_size;
  hash_table_alloc(&db_graph->ht, capacity);
  graph_info_alloc(&db_graph->ginfo);
  return db_graph;
}

void db_graph_dealloc(dBGraph *db_graph)
{
  hash_table_dealloc(&db_graph->ht);
  graph_info_dealloc(&db_graph->ginfo);
}

//
// Add to the de bruijn graph
//

// Note: node may alreay exist in the graph
hkey_t db_graph_find_or_add_node(dBGraph *db_graph, const BinaryKmer bkey,
                                 Colour col)
{
  boolean found;
  hkey_t node = hash_table_find_or_insert(&db_graph->ht, bkey, &found);
  db_node_set_col(db_graph, node, col);
  if(db_graph->covgs != NULL) db_node_increment_coverage(db_graph, node, col);
  return node;
}

// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge(dBGraph *db_graph,
                       hkey_t src_node, hkey_t tgt_node,
                       Orientation src_orient, Orientation tgt_orient)
{
  ConstBinaryKmerPtr src_bkmer = db_node_bkmer(db_graph, src_node);
  ConstBinaryKmerPtr tgt_bkmer = db_node_bkmer(db_graph, tgt_node);

  Nucleotide lhs_nuc, rhs_nuc;
  lhs_nuc = db_node_first_nuc(src_bkmer, src_orient, db_graph->kmer_size);
  rhs_nuc = db_node_last_nuc(tgt_bkmer, tgt_orient, db_graph->kmer_size);

  db_node_set_edge(db_graph, src_node, rhs_nuc, src_orient);
  db_node_set_edge(db_graph, tgt_node, binary_nuc_complement(lhs_nuc),
                   opposite_orientation(tgt_orient));
}

//
// Graph Traversal
//

void db_graph_next_node(const dBGraph *db_graph,
                        const BinaryKmer bkmer, Nucleotide next_nuc,
                        hkey_t *next_node, Orientation *next_orient)
{
  BinaryKmer cpy;
  binary_kmer_assign(cpy, bkmer);
  binary_kmer_left_shift_add(cpy, db_graph->kmer_size, next_nuc);

  BinaryKmer tmp_key;
  db_node_get_key(cpy, db_graph->kmer_size, tmp_key);

  *next_node = hash_table_find(&db_graph->ht, tmp_key);
  *next_orient = db_node_get_orientation(cpy, tmp_key);
}

// Nuc is expected to be already orientated
void db_graph_next_node_orient(const dBGraph *db_graph,
                               const BinaryKmer bkmer, Nucleotide next_nuc,
                               Orientation orient,
                               hkey_t *next_node, Orientation *next_orient)
{
  BinaryKmer tmp_kmer;

  if(orient == forward)
    binary_kmer_assign(tmp_kmer, bkmer);
  else
    binary_kmer_reverse_complement(bkmer, db_graph->kmer_size, tmp_kmer);

  db_graph_next_node(db_graph, tmp_kmer, next_nuc, next_node, next_orient);
}

uint8_t db_graph_next_nodes(const dBGraph *db_graph,
                            const BinaryKmer fw_bkmer, Edges edges,
                            hkey_t nodes[4], BinaryKmer bkmers[4])
{
  // char str[100];
  // binary_kmer_to_str(fw_bkmer, db_graph->kmer_size, str);
  // printf(" :%s [%u]\n", str, edges);

  uint8_t count = 0;
  Edges tmp_edge;
  Nucleotide nuc;
  BinaryKmer bkmer, bkey;

  binary_kmer_assign(bkmer, fw_bkmer);
  binary_kmer_left_shift_one_base(bkmer, db_graph->kmer_size);

  for(tmp_edge = 0x1, nuc = 0; nuc < 4; tmp_edge <<= 1, nuc++)
  {
    if(edges & tmp_edge)
    {
      binary_kmer_set_last_nuc(bkmer, nuc);
      db_node_get_key(bkmer, db_graph->kmer_size, bkey);
      binary_kmer_assign(bkmers[count], bkmer);
      nodes[count] = hash_table_find(&db_graph->ht, bkey);
      count++;

      // binary_kmer_to_str(bkmer, db_graph->kmer_size, str);
      // printf(" ->%s [%i]\n", str, (int)nuc);
    }
  }

  return count;
}

uint8_t db_graph_next_nodes_orient(const dBGraph *db_graph,
                                   const BinaryKmer bkmer, Edges edges,
                                   Orientation orient,
                                   hkey_t nodes[4], BinaryKmer bkmers[4])
{
  BinaryKmer fw_bkmer;
  db_node_oriented_bkmer(bkmer, orient, db_graph->kmer_size, fw_bkmer);
  edges = edges_with_orientation(edges, orient);
  return db_graph_next_nodes(db_graph, fw_bkmer, edges, nodes, bkmers);
}

//
// Pruning
//

static inline void prune_node_without_edges(dBGraph *db_graph, hkey_t node)
{
  if(db_graph->covgs != NULL)
    memset(db_graph->covgs[node], 0, sizeof(db_graph->covgs[node]));

  Colour col;
  for(col = 0; col < NUM_OF_COLOURS; col++)
    db_node_del_col(db_graph, node, col);

  hash_table_delete(&db_graph->ht, node);
  db_graph->ht.unique_kmers--;
}

static void prune_nodes_lacking_flag(hkey_t node, dBGraph *db_graph,
                                     uint64_t *flags)
{
  assert(db_graph->edges != NULL);

  if(bitset_has(flags, node))
  {
    // Check edges
    Orientation orient, next_orient;
    Nucleotide nuc;
    hkey_t next_node;

    Edges keep_edges = db_node_edges(db_graph, node);
    ConstBinaryKmerPtr bkmerptr = db_node_bkmer(db_graph, node);

    for(orient = 0; orient < 2; orient++)
    {
      for(nuc = 0; nuc < 4; nuc++)
      {
        if(edges_has_edge(keep_edges, nuc, orient))
        {
          db_graph_next_node_orient(db_graph, bkmerptr, nuc, orient,
                                    &next_node, &next_orient);
        
          if(!bitset_has(flags, next_node))
          {
            // Next node fails filter - remove edge
            keep_edges = edges_del_edge(keep_edges, nuc, orient);
          }
        }
      }
    }

    db_graph->edges[node] = keep_edges;
  }
  else
  {
    db_node_reset_edges(db_graph, node);
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
  if(db_graph->edges != NULL) {
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
  assert(db_graph->edges != NULL);

  ConstBinaryKmerPtr bkmerptr = db_node_bkmer(db_graph, node);
  hkey_t next_node;
  Orientation or, next_or;
  Nucleotide nuc, lost_nuc;
  Edges remove_edge_mask;

  for(or = 0; or < 2; or++)
  {
    if(edges_with_orientation(edges, or) != 0)
    {
      lost_nuc = db_node_first_nuc(bkmerptr, or, db_graph->kmer_size);

      for(nuc = 0; nuc < 4; nuc++)
      {
        if(edges_has_edge(edges, nuc, or))
        {
          db_graph_next_node_orient(db_graph, bkmerptr, nuc, or,
                                    &next_node, &next_or);

          // Remove edge from next_node to this one
          remove_edge_mask = ~nuc_orient_to_edge(lost_nuc, rev_orient(next_or));

          db_graph->edges[next_node] &= remove_edge_mask;
        }
      }
    }
  }
}

void db_graph_prune_node(dBGraph *db_graph, hkey_t node)
{
  if(db_graph->edges != NULL)
  {
    prune_connected_nodes(db_graph, node, db_graph->edges[node]);
    db_node_reset_edges(db_graph, node);
  }

  prune_node_without_edges(db_graph, node);
}

void db_graph_prune_supernode(dBGraph *db_graph, hkey_t *nodes, size_t len)
{
  if(len == 0) return;
  if(len == 1) db_graph_prune_node(db_graph, nodes[0]);
  else {
    db_graph_prune_node(db_graph, nodes[0]);
    db_graph_prune_node(db_graph, nodes[len-1]);

    size_t i;
    for(i = 1; i+1 < len; i++)
    {
      prune_node_without_edges(db_graph, nodes[i]);
      db_node_reset_edges(db_graph, nodes[i]);
    }
  }
}

//
// Functions applying to whole graph
//

// Remove nodes that are in the hash table but not assigned any colours
static void db_graph_remove_node_if_uncoloured(hkey_t node, dBGraph *db_graph)
{
  Colour c;
  for(c = 0; c < NUM_OF_COLOURS && !db_node_has_col(db_graph, node, c); c++);
  if(c == NUM_OF_COLOURS)
    db_graph_prune_node(db_graph, node);
}

void db_graph_remove_uncoloured_nodes(dBGraph *db_graph)
{
  HASH_TRAVERSE(&db_graph->ht, db_graph_remove_node_if_uncoloured, db_graph);
}

void db_graph_wipe_colour(dBGraph *db_graph, Colour col)
{
  if(db_graph->node_in_cols[0] != NULL)
  {
    size_t size = round_bits_to_words64(db_graph->ht.capacity);
    memset(db_graph->node_in_cols[col], 0, size);
  }

  if(db_graph->covgs != NULL)
  {
    size_t i, capacity = db_graph->ht.capacity;
    Covg (*covgs)[NUM_OF_COLOURS] = db_graph->covgs;
    for(i = 0; i < capacity; i++) covgs[i][col] = 0;
  }

  if(db_graph->col_edges != NULL)
  {
    size_t i, capacity = db_graph->ht.capacity;
    Edges (*col_edges)[NUM_OF_COLOURS] = db_graph->col_edges;
    for(i = 0; i < capacity; i++) col_edges[i][col] = 0;
  }

  db_graph_remove_uncoloured_nodes(db_graph);
}

//
// Kmer paths
//

void db_graph_dump_paths_by_kmer(const dBGraph *db_graph)
{
  const binary_paths_t *paths = &db_graph->pdata;
  path_t path;
  path_alloc(&path);
  uint32_t kmer_size = db_graph->kmer_size;
  char str[100];
  uint64_t node, index;
  Orientation orient;
  boolean printed = false;
  message("\n-------- paths --------\n");

  for(node = 0; node < db_graph->ht.capacity; node++) {
    if(db_graph_node_assigned(db_graph, node)) {
      binary_kmer_to_str(db_node_bkmer(db_graph, node), kmer_size, str);
      for(orient = 0; orient < 2; orient++)
      {
        printed = false;
        index = db_node_paths(db_graph, node, orient);
        while(index != PATH_NULL) {
          if(!printed) { printf("%s:%i\n", str, orient); printed = true; }
          binary_paths_fetch(paths, index, &path);
          binary_paths_dump_path(&path);
          index = path.core.prev;
        }
      }
    }
  }
  path_dealloc(&path);
  message("-----------------------\n\n");
}
