#ifndef DB_NODE_H_
#define DB_NODE_H_

#include <inttypes.h>

#include "graph_typedef.h"
#include "db_graph.h"

//
// Get Binary kmers
//
#define db_node_bkmer(graph,key) ((graph)->ht.table[key])

//
// Get binary kmer key
//
// For a given kmer, get the BinaryKmer 'key':
// the lower of the kmer vs reverse complement of itself
BinaryKmer db_node_get_key(const BinaryKmer kmer, uint32_t kmer_size);

#define db_node_to_str(graph,node,str) \
        binary_kmer_to_str(db_node_bkmer(graph,node), (graph)->kmer_size, (str))

//
// kmer in colours
//
// kset is short for a kmer bitset. A three colour kset is organised like so:
// <0-63:col0><0-63:col1><0-63:col2><64-127:col0><64-127:col1><64-127:col2>
#define ksetw(idx) ((idx)%64)
#define kseto(ncols,col,idx) (((idx)/64) * (ncols) + (col))
#define ksetn(ncols,capacity) (round_bits_to_words(capacity) * (ncols))

#define db_node_has_col(g,hkey,col) \
        bitset2_has((g)->node_in_cols,kseto((g)->num_of_cols,col,hkey),ksetw(hkey))
#define db_node_set_col(g,hkey,col) \
        bitset2_set((g)->node_in_cols,kseto((g)->num_of_cols,col,hkey),ksetw(hkey))
#define db_node_del_col(g,hkey,col) \
        bitset2_del((g)->node_in_cols,kseto((g)->num_of_cols,col,hkey),ksetw(hkey))

//
// Node traversal
//
#define db_node_has_traversed(arr,hkey,or) \
        bitset_has((arr), 2*(hkey)+(or))
#define db_node_set_traversed(arr,hkey,or) \
        bitset_set((arr), 2*(hkey)+(or))

#define db_node_fast_clear_traversed(arr,hkey) \
        bitset_clear_word((arr), 2*(hkey))

//
// Paths
#define db_node_paths(graph,node) ((graph)->kmer_paths[(node)])

//
// Read start (duplicate removal during read loading)
//
#define db_node_has_read_start(graph,hkey,or) \
        bitset_has((graph)->readstrt, 2*(hkey)+(or))
#define db_node_set_read_start(graph,hkey,or) \
        bitset_set((graph)->readstrt, 2*(hkey)+(or))

//
// Orientations
//
#define rev_orient(or) (!(or))
#define opposite_orientation(or) rev_orient(or)

#define db_node_get_orientation(bkmer,bkey) \
        (binary_kmers_are_equal((bkmer), (bkey)) ? FORWARD : REVERSE)

#define db_node_oriented_bkmer(bkmer,or,ksize) \
        (or == FORWARD ? bkmer : binary_kmer_reverse_complement(bkmer,ksize))

#define db_node_first_nuc(bkmer,or,k) \
  ((or) == FORWARD ? binary_kmer_first_nuc((bkmer),(k)) \
      : binary_nuc_complement(binary_kmer_last_nuc(bkmer)))

#define db_node_last_nuc(bkmer,or,k) \
  ((or) == FORWARD ? binary_kmer_last_nuc(bkmer) \
      : binary_nuc_complement(binary_kmer_first_nuc(bkmer,(k))))

//
// Edges
//

// Orientation is 0(FORWARD) or 1(REVERSE), shifted left (<<) 2 gives 0 or 4
#define nuc_orient_to_edge(n,or)    (1 << ((n) + ((or)<<2)))

#define edges_has_edge(edges,n,or)  (((edges) >> ((n) + ((or)<<2))) & 0x1)
#define edges_set_edge(edges,n,or)  ((edges) | nuc_orient_to_edge(n,or))
#define edges_del_edge(edges,n,or)  ((edges) &~nuc_orient_to_edge(n,or))

// shift right by 0 or 4, then AND with 0xf
#define edges_with_orientation(edges,or) (((edges) >> ((or)<<2)) & 0xf)

#define edges_get_outdegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,or))

#define edges_get_indegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,rev_orient(or)))


// #define db_node_has_edge(graph,hkey,nuc,or) \
//         edges_has_edge((graph)->edges[hkey],(nuc),(or))
// #define db_node_del_edge(graph,hkey,nuc,or) \
//         ((graph)->edges[hkey] = edges_del_edge((graph)->edges[hkey],(nuc),(or)))
// #define db_node_set_edge(graph,hkey,nuc,or) \
//         ((graph)->edges[hkey] = edges_set_edge((graph)->edges[hkey],(nuc),(or)))

#define db_node_edges(graph,hkey) db_node_col_edges(graph,0,hkey)
#define db_node_reset_edges(graph,hkey) (db_node_edges(graph,hkey) = 0)

#define db_node_col_edges(graph,col,hkey) \
        ((graph)->col_edges[(hkey)*(graph)->num_edge_cols + (col)])

#define db_node_col_edges_union(graph,hkey) \
        edges_get_union((graph)->col_edges+(hkey)*(graph)->num_edge_cols, \
                        (graph)->num_edge_cols)

#define db_node_set_col_edge(graph,col,hkey,nuc,or) \
        (db_node_col_edges(graph,col,hkey) \
           = edges_set_edge(db_node_col_edges(graph,col,hkey), (nuc), (or)))

Edges edges_get_union(const Edges *edges, size_t num);

boolean edges_has_precisely_one_edge(Edges edges, Orientation orientation,
                                     Nucleotide *nucleotide);

//
// Coverages
//

#define db_node_get_covg(graph,hkey,col) \
        ((graph)->col_covgs[hkey*(graph)->num_of_cols+(col)])

#define db_node_set_covg(graph,hkey,col,covg) \
        (db_node_get_covg(graph,hkey,col) = (covg))

#define db_node_zero_covgs(graph,hkey) \
        memset((graph)->col_covgs + (hkey)*(graph)->num_of_cols, 0, \
               (graph)->num_of_cols * sizeof(Covg))

#define db_node_col_covg(graph,node,colour) \
        ((graph)->col_covgs[(node)*(graph)->num_of_cols + (colour)])

void db_node_add_col_covg(dBGraph *graph, hkey_t hkey, Colour col, long update);
void db_node_increment_coverage(dBGraph *graph, hkey_t hkey, Colour col);

Covg db_node_sum_covg_of_colours(const dBGraph *graph, hkey_t hkey,
                                 Colour first_col, int num_of_cols);

//
// dBNodeBuffer
//
// We might have fewer cache misses if we used this data structure
typedef struct {
  hkey_t node;
  Orientation orient;
} dBNode;

typedef struct
{
  // hkey_t *nodes;
  // Orientation *orients;
  dBNode *data;
  size_t len, capacity;
} dBNodeBuffer;

void db_node_buf_alloc(dBNodeBuffer *buf, size_t capacity);
void db_node_buf_dealloc(dBNodeBuffer *buf);
void db_node_buf_ensure_capacity(dBNodeBuffer *buf, size_t capacity);
void db_node_buf_safe_add(dBNodeBuffer *buf, hkey_t node, Orientation orient);

void db_nodes_to_str(const dBNode *nodes, size_t num,
                     const dBGraph *db_graph, char *str);

#endif /* DB_NODE_H_ */
