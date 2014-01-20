#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <inttypes.h>

#include "string_buffer.h"

#include "cortex_types.h"
#include "hash_table.h"
#include "graph_info.h"
#include "path_store.h"

//
// Graph
//
typedef struct
{
  HashTable ht;
  // num_edge_cols is how many edges are stored per node: 1 or num_of_cols
  const size_t kmer_size;
  const size_t num_of_cols; // How many colours malloc'd for node_in_cols,col_covgs,ginfo
  const size_t num_edge_cols; // How many colours malloc'd for col_edges
  size_t num_of_cols_used; // how many colours currently used

  volatile uint8_t *bktlocks;

  // Array of GraphInfo objects, one per colour
  GraphInfo *ginfo;

  // Optional fields:

  // Colour specific arrays
  // cast to 2d array with:
  // Edges (*col_edges)[graph->num_of_cols]
  //   = (Edges (*)[graph->num_of_cols])graph->col_edges;
  // Covg (*col_covgs)[graph->num_of_cols]
  //   = (Covg (*)[graph->num_of_cols])graph->col_covgs;
  // then access with col_edges[hkey][col]
  Edges *col_edges; // [hkey*num_of_colours + col] or [hkey][col]
  Covg *col_covgs; // [hkey*num_of_colours + col] or [hkey][col]

  // [hkey/64][col] >> hkey%64
  // [num_of_colours*hkey/64+col] >> hkey%64
  uint8_t *node_in_cols;

  // path data
  PathStore pdata;
  // DEV: remove volatile on kmer_paths
  volatile PathIndex *kmer_paths;
  volatile uint8_t *path_kmer_locks;

  // Loading reads
  uint8_t *readstrt;
} dBGraph;

#define db_graph_node_assigned(graph,hkey) HASH_ENTRY_ASSIGNED((graph)->ht.table[hkey])

void db_graph_alloc(dBGraph *db_graph, size_t kmer_size,
                    size_t num_of_cols, size_t num_edge_cols,
                    uint64_t capacity);

void db_graph_realloc(dBGraph *graph, size_t num_of_cols, size_t num_edge_cols);

void db_graph_dealloc(dBGraph *db_graph);

// Get an oriented bkmer
#define db_graph_oriented_bkmer(graph,hkey,or) \
        db_node_oriented_bkmer(db_node_bkmer(graph,hkey),or,(graph)->kmer_size)

//
// Add to the de bruijn graph
//

// Thread safe
// Note: node may alreay exist in the graph
hkey_t db_graph_find_or_add_node_mt(dBGraph *db_graph, BinaryKmer bkey, Colour col);

// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge(dBGraph *db_graph, Colour colour,
                       hkey_t src_node, hkey_t tgt_node,
                       Orientation src_orient, Orientation tgt_orient);

// Thread safe
// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge_mt(dBGraph *db_graph, Colour colour,
                          hkey_t src_node, hkey_t tgt_node,
                          Orientation src_orient, Orientation tgt_orient);

// For debugging + healthcheck
void db_graph_check_edges(const dBGraph *db_graph,
                          hkey_t src_node, hkey_t tgt_node,
                          Orientation src_orient, Orientation tgt_orient);

//
// Graph Traversal
//

void db_graph_next_node(const dBGraph *db_graph, const BinaryKmer node_bkey,
                        Nucleotide next_nuc, Orientation orient,
                        hkey_t *next_node, Orientation *next_orient);

// edges are forward+reverse, db_graph_next_nodes orients them
// fw_nucs is the nuc you would add when walking forward
size_t db_graph_next_nodes(const dBGraph *db_graph, const BinaryKmer node_bkey,
                           Orientation orient, Edges edges,
                           hkey_t nodes[4], Orientation orients[4],
                           Nucleotide fw_nucs[4]);

// Check kmer size of a file
void db_graph_check_kmer_size(size_t kmer_size, const char *path);

//
// Healthcheck
//

void db_graph_healthcheck(const dBGraph *db_graph);

//
// Functions applying to whole graph
//
void db_graph_wipe_colour(dBGraph *db_graph, Colour col);

// Add edges between all kmers with k-1 bases overlapping
void db_graph_add_all_edges(dBGraph *db_graph);

void db_graph_dump_paths_by_kmer(const dBGraph *db_graph);

// Get a random node from the graph
// call seed_random() before any calls to this function please
hkey_t db_graph_rand_node(const dBGraph *db_graph);

//
// Printing
//
void db_graph_print_kmer2(BinaryKmer bkmer, Covg *covgs, Edges *edges,
                          size_t num_of_cols, size_t kmer_size, FILE *fout);

void db_graph_print_kmer(hkey_t node, dBGraph *db_graph, FILE *fout);

#endif /* DB_GRAPH_H_ */
