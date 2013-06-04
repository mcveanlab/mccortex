#include "global.h"
#include "db_graph.h"

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
  if(db_graph->edges != NULL) free(db_graph->edges);
  if(db_graph->status != NULL) free(db_graph->status);
  if(db_graph->bkmer_in_cols != NULL) free(db_graph->bkmer_in_cols);
}

void db_graph_get_next_node_bkmer(const dBGraph *db_graph,
                                  const BinaryKmer bkmer, Nucleotide next_nuc,
                                  hkey_t *next_node, Orientation *next_orient)
{
  BinaryKmer cpy;
  binary_kmer_assign(cpy, bkmer);
  binary_kmer_left_shift_one_base(cpy, db_graph->kmer_size);
  binary_kmer_set_last_nuc(cpy, next_nuc);

  BinaryKmer tmp_key;
  db_node_get_key(cpy, db_graph->kmer_size, tmp_key);

  *next_node = hash_table_find(&db_graph->ht, tmp_key);
  *next_orient = db_node_get_orientation(cpy, tmp_key);
}

void db_graph_get_next_node(const dBGraph *db_graph,
                            const BinaryKmer bkmer, Orientation orient,
                            Nucleotide next_nuc,
                            hkey_t *next_node, Orientation *next_orient)
{
  BinaryKmer tmp_kmer;

  if(orient == forward)
    binary_kmer_assign(tmp_kmer, bkmer);
  else
    binary_kmer_reverse_complement(bkmer, db_graph->kmer_size, tmp_kmer);

  db_graph_get_next_node_bkmer(db_graph, tmp_kmer, next_nuc,
                               next_node, next_orient);
}

int db_graph_get_next_nodes(const dBGraph *db_graph,
                            const BinaryKmer bkmer,
                            Orientation orient, Edges edges,
                            hkey_t nodes[4], BinaryKmer bkmers[4],
                            Orientation orients[4])
{
  int num_of_kmers = db_node_next_bkmers(bkmer, orient, edges,
                                         db_graph->kmer_size, bkmers);

  int i;
  for(i = 0; i < num_of_kmers; i++)
  {
    BinaryKmer tmp_key;
    db_node_get_key(bkmers[i], db_graph->kmer_size, tmp_key);

    nodes[i] = hash_table_find(&db_graph->ht, tmp_key);
    orients[i] = db_node_get_orientation(bkmers[i], tmp_key);
  }

  return num_of_kmers;
}

// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge(dBGraph *db_graph,
                       hkey_t src_node, hkey_t tgt_node,
                       Orientation src_orient, Orientation tgt_orient)
{
  BinaryKmerPtr src_bkmer = db_graph_bkmer(db_graph, src_node);
  BinaryKmerPtr tgt_bkmer = db_graph_bkmer(db_graph, tgt_node);

  Nucleotide lhs_nuc, rhs_nuc;
  lhs_nuc = db_node_first_nuc(src_bkmer, src_orient, db_graph->kmer_size);
  rhs_nuc = db_node_last_nuc(tgt_bkmer, tgt_orient, db_graph->kmer_size);

  db_node_set_edge(db_graph, src_node, rhs_nuc, src_orient);
  db_node_set_edge(db_graph, tgt_node, binary_nucleotide_complement(lhs_nuc),
                   opposite_orientation(tgt_orient));
}

//
// Pruning
//

static inline void prune_node_without_edges(dBGraph *db_graph, hkey_t node)
{
  if(db_graph->covgs != NULL)
    memset(db_graph->covgs[node], 0, sizeof(db_graph->covgs[node]));

  if(db_graph->status != NULL)
    db_graph->status[node] = EFLAG_ZERO;

  Colour col;
  for(col = 0; col < NUM_OF_COLOURS; col++)
    db_graph_bkmer_del_col(db_graph, node, col);

  hash_table_delete(&db_graph->ht, node);
  db_graph->ht.unique_kmers--;
}

static void prune_nodes_lacking_flag(hkey_t node, dBGraph *db_graph, uint8_t flag)
{
  assert(db_graph->edges != NULL);

  if(db_node_has_flag(db_graph, node, flag))
  {
    // Check edges
    Orientation orient, next_orient;
    Nucleotide nuc;
    hkey_t next_node;

    Edges keep_edges = db_node_get_edges(db_graph, node);
    BinaryKmerPtr bkmer = db_graph_bkmer(db_graph, node);

    for(orient = 0; orient < 2; orient++)
    {
      for(nuc = 0; nuc < 4; nuc++)
      {
        if(edges_has_edge(keep_edges, nuc, orient))
        {
          db_graph_get_next_node(db_graph, bkmer, orient, nuc,
                                 &next_node, &next_orient);
        
          if(!db_node_has_flag(db_graph, next_node, flag))
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

static void prune_nodes_lacking_flag2(hkey_t node, dBGraph *db_graph, uint8_t flag)
{
  if(!db_node_has_flag(db_graph, node, flag))
    prune_node_without_edges(db_graph, node);
}

// If element_has_flag(node, flag) is not true, reset flags
// Remove edges to nodes where !db_node_has_flag(node, flag)
void db_graph_prune_nodes_lacking_flag(dBGraph *db_graph, uint8_t flag)
{
  if(db_graph->edges != NULL) {
    HASH_TRAVERSE(&db_graph->ht, prune_nodes_lacking_flag, db_graph, flag);
  }
  else {
    HASH_TRAVERSE(&db_graph->ht, prune_nodes_lacking_flag2, db_graph, flag);
  }
}

// Removed edges from nodes that connected to the given `node`
static void prune_connected_nodes(dBGraph *db_graph, hkey_t node,
                                  Orientation orient, Edges edges)
{
  assert(db_graph->edges != NULL);

  hkey_t next_node;
  Orientation next_orient;
  Nucleotide nuc, lost_nuc;
  Edges remove_edge_mask;

  if(db_node_edges_in_orientation(edges, orient) == 0) return;

  lost_nuc = db_node_first_nuc(db_graph_bkmer(db_graph, node),
                               orient, db_graph->kmer_size);

  BinaryKmerPtr bkmer = db_graph_bkmer(db_graph, node);

  for(nuc = 0; nuc < 4; nuc++)
  {
    if(edges_has_edge(edges, nuc, orient))
    {

      db_graph_get_next_node(db_graph, bkmer, orient, nuc,
                             &next_node, &next_orient);

      // Remove edge from next_node to this one
      remove_edge_mask = ~nuc_orient_to_edge(lost_nuc,
                                             opposite_orientation(next_orient));

      db_graph->edges[next_node] &= remove_edge_mask;
    }
  }
}

void db_graph_prune_node(dBGraph *db_graph, hkey_t node)
{
  if(db_graph->edges != NULL)
  {
    Edges edges = db_graph->edges[node];
    Orientation orient;
    for(orient = 0; orient < 2; orient++) {
      prune_connected_nodes(db_graph, node, orient, edges);
    }
    db_node_reset_edges(db_graph, node);
  }

  prune_node_without_edges(db_graph, node);
}

// Edges are required
void db_graph_prune_nodes(dBGraph *db_graph, hkey_t *nodes, size_t len,
                          boolean is_supernode)
{
  if(db_graph->edges == NULL)
  {
    size_t i;
    for(i = 0; i < len; i++)
      prune_node_without_edges(db_graph, nodes[i]);
    return;
  }

  // If nodes have in-degree==1 && out-degree==1 and they're in the middle of
  // the contig we don't need to look up their neighbours in the graph
  // in-degree==1 && out-degree==1 is guaranteed for supernodes

  db_graph_prune_node(db_graph, nodes[0]);
  db_graph_prune_node(db_graph, nodes[len-1]);

  size_t i;

  if(!is_supernode)
  {
    for(i = 1; i+1 < len; i++)
    {
      Edges edges = db_node_get_edges(db_graph, nodes[i]);

      if(edges_get_outdegree(edges, forward) > 1) {
        prune_connected_nodes(db_graph, nodes[i], forward, edges);
      }
      if(edges_get_outdegree(edges, reverse) > 1) {
        prune_connected_nodes(db_graph, nodes[i], reverse, edges);
      }
    }
  }

  for(i = 1; i+1 < len; i++)
    prune_node_without_edges(db_graph, nodes[i]);
}

//
// Functions applying to whole graph
//

// Remove nodes that are in the hash table but not assigned any colours
static void db_graph_remove_node_if_uncoloured(hkey_t node, dBGraph *db_graph)
{
  Colour col = 0;
  while(col < NUM_OF_COLOURS && !db_graph_bkmer_has_col(db_graph, node, col)) {
    col++;
  }
  if(col == NUM_OF_COLOURS) db_graph_prune_node(db_graph, node);
}

void db_graph_remove_uncoloured_nodes(dBGraph *db_graph)
{
  HASH_TRAVERSE(&db_graph->ht, db_graph_remove_node_if_uncoloured, db_graph);
}

void db_graph_wipe_colour(dBGraph *db_graph, Colour col)
{
  size_t bitfield_size = db_graph_sizeof_bkmer_bitset(db_graph);
  memset(db_graph->bkmer_in_cols[col], 0, bitfield_size);

  if(db_graph->covgs != NULL)
  {
    size_t i, capacity = db_graph->ht.capacity;
    Covg (*covgs)[NUM_OF_COLOURS] = db_graph->covgs;
    for(i = 0; i < capacity; i++) covgs[i][col] = 0;
  }

  db_graph_remove_uncoloured_nodes(db_graph);
}
