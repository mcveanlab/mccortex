#ifndef DB_NODE_H_
#define DB_NODE_H_

#include <inttypes.h>
#include "bit_macros.h"

#include "cortex_types.h"
#include "db_graph.h"

#define db_nodes_match(n1,n2) ((n1).key == (n2).key && (n1).orient == (n2).orient)

//
// Get Binary kmers
//
#define db_node_bkmer(graph,key) ((graph)->ht.table[key])

static inline BinaryKmer db_node_get_bkmer(const dBGraph *db_graph, hkey_t hkey) {
  return db_graph->ht.table[hkey];
}

//
// Get binary kmer key
//
// For a given kmer, get the BinaryKmer 'key':
// the lower of the kmer vs reverse complement of itself
BinaryKmer db_node_get_key(const BinaryKmer kmer, size_t kmer_size);

#define db_node_to_str(graph,node,str) \
        binary_kmer_to_str(db_node_bkmer(graph,node), (graph)->kmer_size, (str))

//
// kmer in colours
//

// kset is short for a kmer bitset. A three colour kset is organised like so:
// <0-63:col0><0-63:col1><0-63:col2><64-127:col0><64-127:col1><64-127:col2>

/* Offset in word */
#define kseto(arr,idx) ((idx)%(sizeof(*arr)*8))
/* word index */
#define ksetw(arr,ncols,col,idx) (((idx)/(sizeof(*arr)*8))*(ncols)+(col))

static inline boolean db_node_has_col(const dBGraph *graph, hkey_t hkey, size_t col) {
  return bitset2_get(graph->node_in_cols,
                     ksetw(graph->node_in_cols,graph->num_of_cols,col,hkey),
                     kseto(graph->node_in_cols,hkey));
}

static inline void db_node_set_col(const dBGraph *graph, hkey_t hkey, size_t col) {
  bitset2_set(graph->node_in_cols,
              ksetw(graph->node_in_cols,graph->num_of_cols,col,hkey),
              kseto(graph->node_in_cols,hkey));
}

static inline void db_node_del_col(const dBGraph *graph, hkey_t hkey, size_t col) {
  bitset2_del(graph->node_in_cols,
              ksetw(graph->node_in_cols,graph->num_of_cols,col,hkey),
              kseto(graph->node_in_cols,hkey));
}

static inline void db_node_cpy_col(const dBGraph *graph, hkey_t hkey,
                                   size_t col, uint8_t bit) {
  bitset2_cpy(graph->node_in_cols,
              ksetw(graph->node_in_cols,graph->num_of_cols,col,hkey),
              kseto(graph->node_in_cols,hkey),bit);
}

// Threadsafe
static inline void db_node_set_col_mt(const dBGraph *graph,
                                      hkey_t hkey, size_t col) {
  bitset2_set_mt(graph->node_in_cols,
                 ksetw(graph->node_in_cols,graph->num_of_cols,col,hkey),
                 kseto(graph->node_in_cols,hkey));
}


//
// Node traversal
//
#define db_node_has_traversed(arr,node) \
        bitset_get((arr), 2*((node).key)+((node).orient))
#define db_node_set_traversed(arr,node) \
        bitset_set((arr), 2*((node).key)+((node).orient))

#define db_node_fast_clear_traversed(arr,hkey) \
        bitset_clear_word((arr), 2*(hkey))

//
// Paths
#define db_node_paths(graph,node) ((graph)->kmer_paths[(node)])

//
// Orientations
//
#define rev_orient(or) (!(or))
#define opposite_orientation(or) rev_orient(or)

#define db_node_get_orientation(bkmer,bkey) \
        (binary_kmers_are_equal((bkmer), (bkey)) ? FORWARD : REVERSE)

#define db_node_oriented_bkmer(bkmer,or,ksize) \
        (or == FORWARD ? bkmer : binary_kmer_reverse_complement(bkmer,ksize))

#define db_node_first_nuc(bkmer,or,ksize) \
  ((or) == FORWARD ? binary_kmer_first_nuc(bkmer,ksize) \
                   : dna_nuc_complement(binary_kmer_last_nuc(bkmer)))

#define db_node_last_nuc(bkmer,or,ksize) \
  ((or) == FORWARD ? binary_kmer_last_nuc(bkmer) \
                   : dna_nuc_complement(binary_kmer_first_nuc(bkmer,ksize)))

static inline dBNode db_node_reverse(dBNode node) {
  node.orient = !node.orient;
  return node;
}

//
// Edges
//

// Orientation is 0(FORWARD) or 1(REVERSE), shifted left (<<) 2 gives 0 or 4
#define nuc_orient_to_edge(n,or)    ((Edges)(1 << ((n) + ((or)<<2))))

#define edges_has_edge(edges,n,or)  (((edges) >> ((n) + ((or)<<2))) & 0x1)
#define edges_set_edge(edges,n,or)  ((edges) | nuc_orient_to_edge(n,or))
#define edges_del_edge(edges,n,or)  ((edges) &~nuc_orient_to_edge(n,or))

// shift right by 0 or 4, then AND with 0xf
#define edges_with_orientation(edges,or) (((edges) >> ((or)<<2)) & 0xf)

#define edges_get_outdegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,or))

#define edges_get_indegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,rev_orient(or)))

Edges edges_get_union(const Edges *edges, size_t num);

boolean edges_has_precisely_one_edge(Edges edges, Orientation orientation,
                                     Nucleotide *nucleotide);

//
// dBNode Edges
//

#define db_node_edges(graph,col,hkey) \
        ((graph)->col_edges[(hkey)*(graph)->num_edge_cols + (col)])

static inline Edges db_node_get_edges(const dBGraph *graph, Colour col, hkey_t hkey) {
  return db_node_edges(graph, col, hkey);
}

static inline Edges db_node_get_edges_union(const dBGraph *graph, hkey_t hkey) {
  return edges_get_union(graph->col_edges + hkey * graph->num_edge_cols,
                         graph->num_edge_cols);
}

Edges db_node_oriented_edges_in_col(dBNode node, size_t col,
                                    const dBGraph *db_graph);

#define db_node_zero_edges(graph,hkey) \
        memset((graph)->col_edges + (hkey)*(graph)->num_edge_cols, 0, \
               (graph)->num_edge_cols * sizeof(Edges))

#define db_node_set_col_edge(graph,col,hkey,nuc,or) \
        (db_node_edges(graph,col,hkey) \
           = edges_set_edge(db_node_get_edges(graph,col,hkey),nuc,or))

#define db_node_set_col_edge_mt(graph,col,hkey,nuc,or) \
        __sync_or_and_fetch(&db_node_edges(graph,col,hkey), nuc_orient_to_edge(nuc,or))

// kmer_col_edge_str should be 9 chars long
// Return pointer to kmer_col_edge_str
char* db_node_get_edges_str(Edges edges, char* kmer_col_edge_str);

//
// Coverages
//

#define db_node_covg(graph,hkey,col) \
        ((graph)->col_covgs[(hkey)*(graph)->num_of_cols+(col)])

#define db_node_zero_covgs(graph,hkey) \
        memset((graph)->col_covgs + (hkey)*(graph)->num_of_cols, 0, \
               (graph)->num_of_cols * sizeof(Covg))

#define db_node_col_covg(graph,colour,node) \
        ((graph)->col_covgs[(node)*(graph)->num_of_cols + (colour)])

void db_node_add_col_covg(dBGraph *graph, hkey_t hkey, Colour col, Covg update);
void db_node_increment_coverage(dBGraph *graph, hkey_t hkey, Colour col);

// Thread safe, overflow safe, coverage increment
void db_node_increment_coverage_mt(dBGraph *graph, hkey_t hkey, Colour col);

Covg db_node_sum_covg(const dBGraph *graph, hkey_t hkey);

//
// dBNodeBuffer
//

// The following is defined by the create_objbuf MACRO
// typedef struct {
//   dBNode *data;
//   size_t len, capacity;
// } dBNodeBuffer;

// void db_node_buf_alloc(dBNodeBuffer *buf, size_t capacity);
// void db_node_buf_dealloc(dBNodeBuffer *buf);
// void db_node_buf_ensure_capacity(dBNodeBuffer *buf, size_t capacity);
// void db_node_buf_add(dBNodeBuffer *buf, dBNode node);

#include "objbuf_macro.h"
create_objbuf(db_node_buf,dBNodeBuffer,dBNode)

#define db_node_buf_safe_add(buf,node,or) {\
  dBNode n = {.key=node,.orient=or};  db_node_buf_add(buf,n); }

//
// Print array of dBNode
//

void db_nodes_to_str(const dBNode *nodes, size_t num,
                     const dBGraph *db_graph, char *str);

void db_nodes_print(const dBNode *nodes, size_t num,
                    const dBGraph *db_graph, FILE *out);

void db_nodes_gzprint(const dBNode *nodes, size_t num,
                      const dBGraph *db_graph, gzFile out);

#endif /* DB_NODE_H_ */
