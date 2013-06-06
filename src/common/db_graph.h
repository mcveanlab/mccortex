#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <inttypes.h>

#include "binary_kmer.h"
#include "hash_table.h"

typedef uint8_t Colour;
typedef uint8_t Edges;
typedef uint32_t Covg;

#include "graph_info.h"
#include "binary_paths.h"

typedef enum
{
  forward = 0,
  reverse = 1
} Orientation;

typedef BinaryKmerPtr Key;

typedef struct
{
  HashTable ht;
  uint32_t kmer_size;
  // Optional fields:
  Covg (*covgs)[NUM_OF_COLOURS]; // [hkey][col]
  Edges *edges;
  Edges (*col_edges)[NUM_OF_COLOURS]; // [hkey][col]
  uint64_t *bkmer_in_cols[NUM_OF_COLOURS];
  // path data
  uint64_t *kmer_paths;
  binary_paths_t pdata;
  // Visited
  uint64_t *visited;
  // Loading reads
  uint64_t *readstrt;
  // Info stored here:
  GraphInfo ginfo;
} dBGraph;

// Basic operations on the graph nodes
#include "db_node.h"

dBGraph* db_graph_alloc(dBGraph *db_graph, uint32_t kmer_size, uint64_t capacity);
void db_graph_dealloc(dBGraph *db_graph);

// Get an oriented bkmer
#define db_graph_oriented_bkmer(graph,hkey,or,result) \
        db_node_oriented_bkmer(db_node_bkmer(graph,hkey), or, \
                               (graph)->kmer_size, result)

//
// Add to the de bruijn graph
//

// Note: node may alreay exist in the graph
hkey_t db_graph_find_or_add_node(dBGraph *db_graph, const BinaryKmer bkey,
                                 Colour col);

// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge(dBGraph *db_graph,
                       hkey_t src_node, hkey_t tgt_node,
                       Orientation src_orient, Orientation tgt_orient);

//
// Graph Traversal
//

void db_graph_next_node(const dBGraph *db_graph,
                        const BinaryKmer bkmer, Nucleotide next_nuc,
                        hkey_t *next_node, Orientation *next_orient);

void db_graph_next_node_orient(const dBGraph *db_graph,
                               const BinaryKmer bkmer, Nucleotide next_nuc,
                               Orientation orient,
                               hkey_t *next_node, Orientation *next_orient);

uint8_t db_graph_next_nodes(const dBGraph *db_graph,
                            const BinaryKmer fw_bkmer, Edges edges,
                            hkey_t nodes[4], BinaryKmer bkmers[4]);

uint8_t db_graph_next_nodes_orient(const dBGraph *db_graph,
                                   const BinaryKmer bkmer, Edges edges,
                                   Orientation orient,
                                   hkey_t nodes[4], BinaryKmer bkmers[4]);

//
// Pruning
//

void db_graph_prune_nodes_lacking_flag(dBGraph *graph, uint64_t *flags);

void db_graph_prune_node(dBGraph *db_graph, hkey_t node);
void db_graph_prune_supernode(dBGraph *db_graph, hkey_t *nodes, size_t len);

//
// Functions applying to whole graph
//
void db_graph_remove_uncoloured_nodes(dBGraph *db_graph);

void db_graph_wipe_colour(dBGraph *db_graph, Colour col);

#endif /* DB_GRAPH_H_ */
