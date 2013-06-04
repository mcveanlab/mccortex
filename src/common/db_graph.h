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
  uint8_t *status;
  uint64_t *bkmer_in_cols[NUM_OF_COLOURS];
  // path data
  uint64_t *kmer_paths;
  binary_paths_t pdata;
  // Traversal
  uint64_t *traverse_fw, *traverse_rv;
  // Loading nodes
  uint64_t *readstrt_fw, *readstrt_rv;
  // Info stored here:
  GraphInfo ginfo;
} dBGraph;

// Basic operations on the graph nodes
#include "db_node.h"

dBGraph* db_graph_alloc(dBGraph *db_graph, uint32_t kmer_size, uint64_t capacity);
void db_graph_dealloc(dBGraph *db_graph);

#define db_graph_bkmer(graph,key) ((const BinaryKmerPtr)((graph)->ht.table[key]))

//
// kmer in colours
//
#define db_graph_bkmer_has_col(graph,hkey,col) \
        bitset_has((graph)->bkmer_in_cols[col], hkey)
#define db_graph_bkmer_set_col(graph,hkey,col) \
        bitset_set((graph)->bkmer_in_cols[col], hkey)
#define db_graph_bkmer_del_col(graph,hkey,col) \
        bitset_del((graph)->bkmer_in_cols[col], hkey)

//
// Node traversal
//
#define db_node_has_traversed(graph,hkey,or) \
        bitset_has((or)==forward?(graph)->traverse_fw:(graph)->traverse_rv,hkey)
#define db_node_set_traversed(graph,hkey,or) \
        bitset_set((or)==forward?(graph)->traverse_fw:(graph)->traverse_rv,hkey)

#define db_node_fast_clear_traversed(graph,hkey) \
        {bitset_clear_word((graph)->traverse_fw,hkey); \
         bitset_clear_word((graph)->traverse_fw,hkey); }

//
// Paths
#define db_graph_kmer_path(graph,node,or) ((graph)->kmer_paths[2*node+or])

//
// Read start (duplicate removal during read loading)
//
#define db_node_has_read_start(graph,hkey,or) \
        bitset_has((or)==forward?(graph)->readstrt_fw:(graph)->readstrt_rv,hkey)
#define db_node_set_read_start_status(graph,hkey,or) \
        bitset_set((or)==forward?(graph)->readstrt_fw:(graph)->readstrt_rv,hkey)

// Get number 
#define db_graph_sizeof_bkmer_bitset(graph) \
        round_bits_to_words64((graph)->ht.capacity)

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
                            hkey_t nodes[4], BinaryKmer bkmers[4],
                            Orientation orients[4]);

uint8_t db_graph_next_nodes_orient(const dBGraph *db_graph,
                                   const BinaryKmer bkmer, Edges edges,
                                   Orientation orient,
                                   hkey_t nodes[4], BinaryKmer bkmers[4],
                                   Orientation orients[4]);

// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge(dBGraph *db_graph,
                       hkey_t src_node, hkey_t tgt_node,
                       Orientation src_orient, Orientation tgt_orient);

//
// Pruning
//

void db_graph_prune_nodes_lacking_flag(dBGraph *graph, uint8_t flag);

void db_graph_prune_node(dBGraph *db_graph, hkey_t node);
void db_graph_prune_nodes(dBGraph *db_graph, hkey_t *nodes, size_t len,
                          boolean is_supernode);

void db_graph_prune_node_in_colour(dBGraph *db_graph, hkey_t node, Colour col);

//
// Functions applying to whole graph
//
void db_graph_remove_uncoloured_nodes(dBGraph *db_graph);

void db_graph_wipe_colour(dBGraph *db_graph, Colour col);

#define db_graph_set_all_node_statuses(graph,status) \
        memset((graph)->status, status, (graph)->capacity);

#endif /* DB_GRAPH_H_ */
