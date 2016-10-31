#ifndef DB_NODE_H_
#define DB_NODE_H_

#include <inttypes.h>
#include "bit_array/bit_macros.h"
#include "htslib/khash.h"

#include "cortex_types.h"
#include "db_graph.h"
#include "util.h"

#define DB_NODE_INIT {.key = HASH_NOT_FOUND, .orient = FORWARD}

#define db_nodes_are_equal(n1,n2) ((n1).key == (n2).key && (n1).orient == (n2).orient)

static inline uint64_t db_node_hash(dBNode node) {
  uint64_t d = (node.key << 1) | node.orient;
  return kh_int64_hash_func(d);
}

//
// Get Binary kmers
//
static inline BinaryKmer db_node_get_bkey(const dBGraph *db_graph, hkey_t hkey) {
  return hash_table_fetch(&db_graph->ht, hkey);
}

// Get an oriented bkmer
#define db_node_oriented_bkmer(graph,node) \
        bkmer_oriented_bkmer(db_node_get_bkey(graph,(node).key), \
                             (node).orient, (graph)->kmer_size)

//
// kmer in colours
//

// kset is short for a kmer bitset. A three colour kset is organised like so:
// <0-63:col0><0-63:col1><0-63:col2><64-127:col0><64-127:col1><64-127:col2>
// This is so we can wipe a single colour by zero-ing whole bytes, whilst still
// keeping the same node in multiple colours close together.

/* Offset in word */
#define kseto(arr,hkey)           ((hkey)%(sizeof(*arr)*8))
/* word index */
#define ksetw(arr,ncols,hkey,col) (((hkey)/(sizeof(*arr)*8))*(ncols)+(col))

static inline bool db_node_in_col(const dBGraph *graph, hkey_t hkey, size_t col)
{
  return graph->node_in_cols == NULL ||
         bitset2_get(graph->node_in_cols,
                     ksetw(graph->node_in_cols,graph->num_of_cols,hkey,col),
                     kseto(graph->node_in_cols,hkey));
}

static inline bool db_node_has_col(const dBGraph *graph, hkey_t hkey, size_t col)
{
  return bitset2_get(graph->node_in_cols,
                     ksetw(graph->node_in_cols,graph->num_of_cols,hkey,col),
                     kseto(graph->node_in_cols,hkey));
}

static inline void db_node_set_col(const dBGraph *graph, hkey_t hkey, size_t col)
{
  bitset2_set(graph->node_in_cols,
              ksetw(graph->node_in_cols,graph->num_of_cols,hkey,col),
              kseto(graph->node_in_cols,hkey));
}

static inline void db_node_del_col_mt(const dBGraph *graph, hkey_t hkey, size_t col)
{
  (void)bitset2_del_mt(graph->node_in_cols,
                       ksetw(graph->node_in_cols,graph->num_of_cols,hkey,col),
                       kseto(graph->node_in_cols,hkey));
}

static inline void db_node_or_col(const dBGraph *graph, hkey_t hkey,
                                  size_t col, uint8_t bit)
{
  bitset2_or(graph->node_in_cols,
             ksetw(graph->node_in_cols,graph->num_of_cols,hkey,col),
             kseto(graph->node_in_cols,hkey),bit);
}

// Threadsafe
static inline void db_node_set_col_mt(const dBGraph *graph,
                                      hkey_t hkey, size_t col)
{
  (void)bitset2_set_mt(graph->node_in_cols,
                       ksetw(graph->node_in_cols,graph->num_of_cols,hkey,col),
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
// BinaryKmer + Orientations (bkmer_*)
//

#define bkmer_get_orientation(bkmer,bkey) \
        (binary_kmer_eq((bkmer), (bkey)) ? FORWARD : REVERSE)

#define bkmer_oriented_bkmer(bkmer,or,ksize) \
        (or == FORWARD ? bkmer : binary_kmer_reverse_complement(bkmer,ksize))

#define bkmer_get_first_nuc(bkmer,or,ksize) \
  ((or) == FORWARD ? binary_kmer_first_nuc(bkmer,ksize) \
                   : dna_nuc_complement(binary_kmer_last_nuc(bkmer)))

#define bkmer_get_last_nuc(bkmer,or,ksize) \
  ((or) == FORWARD ? binary_kmer_last_nuc(bkmer) \
                   : dna_nuc_complement(binary_kmer_first_nuc(bkmer,ksize)))

#define bkmer_shift_add_first_nuc(bkey,or,ksize,nuc) \
  ((or) == FORWARD ? binary_kmer_right_shift_add(bkey,ksize,nuc) \
                   : binary_kmer_left_shift_add(bkey,ksize,dna_nuc_complement(nuc)))

#define bkmer_shift_add_last_nuc(bkey,or,ksize,nuc) \
  ((or) == FORWARD ? binary_kmer_left_shift_add(bkey,ksize, nuc) \
                   : binary_kmer_right_shift_add(bkey,ksize,dna_nuc_complement(nuc)))

//
// Orientations
//
#define rev_orient(or) (!(or))
#define opposite_orientation(or) rev_orient(or)

#define db_node_get_first_nuc(node,graph) \
        bkmer_get_first_nuc(db_node_get_bkey(graph,(node).key), (node).orient,\
                            (graph)->kmer_size)

#define db_node_get_last_nuc(node,graph) \
        bkmer_get_last_nuc(db_node_get_bkey(graph,(node).key), (node).orient,\
                           (graph)->kmer_size)

static inline dBNode db_node_reverse(dBNode node) {
  node.orient = !node.orient;
  return node;
}

//
// Edges
//
/*

Edges are coded as:

  v preceeding base on the reverse strand
  acgtACGT
      ^ next base on forward strand

So all on forward strand the coding is:

  TGCAACGT

*/

// Get egdes on forward strand
static inline Edges edges_as_nibble(Edges edges, Orientation orient) {
  if(orient == REVERSE) {
    edges >>= 4;
    edges = rev_nibble_lookup(edges);
  }
  else {
    edges &= 0xf;
  }
  return edges;
}

// Orientation is 0(FORWARD) or 1(REVERSE), shifted left (<<) 2 gives 0 or 4
#define nuc_orient_to_edge(n,or)    ((Edges)(1 << ((n) + ((or)<<2))))

#define edges_has_edge(edges,n,or)  (((edges) >> ((n) + ((or)<<2))) & 0x1)
#define edges_set_edge(edges,n,or)  ((edges) | nuc_orient_to_edge(n,or))
#define edges_del_edge(edges,n,or)  ((edges) &~nuc_orient_to_edge(n,or))

// shift right by 0 or 4, then AND with 0xf
#define edges_with_orientation(edges,or) (((edges) >> ((or)<<2)) & 0xf)
#define edges_mask_orientation(edges,or) ((edges) & (0xf << ((or)<<2)))

#define edges_get_outdegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,or))

#define edges_get_indegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,rev_orient(or)))

Edges edges_get_union(const Edges *edges, size_t num);

bool edges_has_precisely_one_edge(Edges edges, Orientation orientation,
                                  Nucleotide *nucleotide);

// Get edges in hex coding, two characters [0-9a-f] per edge
// 1=>A, 2=>C, 4=>G, 8=>T
// "3b" => [AC] AACTA [ACT]
// Null terminates string
static inline char* edges_to_char(Edges e, char str[3])
{
  static const char digits[16] = "0123456789abcdef";
  str[0] = digits[edges_as_nibble(e, REVERSE)];
  str[1] = digits[edges_as_nibble(e, FORWARD)];
  str[2] = '\0';
  return str;
}

static inline void edges_print(FILE *fout, Edges e)
{
  char estr[3];
  edges_to_char(e, estr);
  fputc(estr[0], fout);
  fputc(estr[1], fout);
}

// str[] must be at least 5 bytes long
static inline char *edges_get_str(Edges fw_edges, char *str)
{
  char bases[] = "ACGT", *start = str;
  int i;
  Edges mask;
  for(i = 0, mask = 1; i < 4; i++, mask<<=1)
    if(fw_edges & mask) *(str++) = bases[i];
  *str = '\0';
  return start;
}



//
// dBNode Edges
//

#define db_node_edges(graph,hkey,col) \
        ((graph)->col_edges[(hkey)*(graph)->num_edge_cols + (col)])

static inline Edges db_node_get_edges(const dBGraph *graph, hkey_t hkey, Colour col) {
  return db_node_edges(graph, hkey, col);
}

static inline Edges db_node_get_edges_union(const dBGraph *graph, hkey_t hkey) {
  return edges_get_union(graph->col_edges + hkey * graph->num_edge_cols,
                         graph->num_edge_cols);
}

// Edges restricted to this colour, only in one direction (node.orient)
Edges db_node_edges_in_col(dBNode node, size_t col, const dBGraph *db_graph);

// Edges in both directions
Edges db_node_both_edges_in_col(hkey_t hkey, size_t col, const dBGraph *db_graph);

// Outdegree of edges restricted to a given colour
#define db_node_outdegree_in_col(node,col,graph) \
        edges_get_outdegree(db_node_edges_in_col(node,col,graph), (node).orient)

#define db_node_indegree_in_col(node,col,graph) \
        db_node_outdegree_in_col(db_node_reverse(node),col,graph)

#define db_node_zero_edges(graph,hkey) \
        memset((graph)->col_edges + (hkey)*(graph)->num_edge_cols, 0, \
               (graph)->num_edge_cols * sizeof(Edges))

#define db_node_set_col_edge(graph,hkey,col,nuc,or) \
        (db_node_edges(graph,hkey,col) \
           = edges_set_edge(db_node_get_edges(graph,hkey,col),nuc,or))

#define db_node_set_col_edge_mt(graph,hkey,col,nuc,or) \
        __sync_or_and_fetch(&db_node_edges(graph,hkey,col), nuc_orient_to_edge(nuc,or))

// kmer_col_edge_str should be 9 chars long
// Return pointer to kmer_col_edge_str
char* db_node_get_edges_str(Edges edges, char *kmer_col_edge_str);

//
// Coverages
//

#define db_node_covg(graph,hkey,col) \
        ((graph)->col_covgs[(hkey)*(graph)->num_of_cols+(col)])

static inline Covg db_node_get_covg(const dBGraph *db_graph,
                                    hkey_t hkey, Colour col) {
  return db_node_covg(db_graph, hkey, col);
}

#define db_node_zero_covgs(graph,hkey) \
        memset((graph)->col_covgs + (hkey)*(graph)->num_of_cols, 0, \
               (graph)->num_of_cols * sizeof(Covg))

void db_node_add_col_covg(dBGraph *graph, hkey_t hkey, Colour col, Covg update);
void db_node_increment_coverage(dBGraph *graph, hkey_t hkey, Colour col);

// Thread safe, overflow safe, coverage increment
void db_node_increment_coverage_mt(dBGraph *graph, hkey_t hkey, Colour col);

static inline Covg db_node_sum_covg(const dBGraph *graph, hkey_t hkey)
{
  const Covg *covgs = &db_node_covg(graph,hkey,0);
  Covg sum_covg = covgs[0];
  size_t c, ncols = graph->num_of_cols;
  for(c = 1; c < ncols; c++) SAFE_SUM_COVG(sum_covg, covgs[c]);
  return sum_covg;
}

//
// dBNodeBuffer
//

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(db_node_buf,dBNodeBuffer,dBNode);

static inline dBNode db_nodes_get(const dBNode *nodes, size_t n, bool fw, size_t i)
{
  return fw ? nodes[i] : db_node_reverse(nodes[n-1-i]);
}

//
// dBNode reversal and shifting
//

// Reverse ordering without changing node orientation
void db_nodes_reverse(dBNode *nlist, size_t n);

// Reverse order and orientations of nodes
void db_nodes_reverse_complement(dBNode *nlist, size_t len);

void db_nodes_left_shift(dBNode *nlist, size_t n, size_t shift);

//
// Print array of dBNode
//

// Get bkey:orient string representation e.g. "AGAGTTTTATC:1".
//   :0 means forward, :1 means reverse
//   `str` must be at least kmer_size+3 chars long
// Returns length in bytes. Null terminates `str`.
size_t db_node_to_str(const dBGraph *db_graph, dBNode node, char *str);

// Null-terminates string
// Returns number of bytes added not including \0
size_t db_nodes_to_str(const dBNode *nodes, size_t num,
                       const dBGraph *db_graph, char *str);

void db_nodes_print(const dBNode *nodes, size_t num,
                    const dBGraph *db_graph, FILE *out);

void db_nodes_gzprint(const dBNode *nodes, size_t num,
                      const dBGraph *db_graph, gzFile out);

// Do not print first k-1 bases => 3 nodes gives 3bp instead of 3+k-1
void db_nodes_gzprint_cont(const dBNode *nodes, size_t num,
                           const dBGraph *db_graph, gzFile out);

// Print:
// 0: AAACCCAAATGCAAACCCAAATGCAAACCCA:1 TGGGTTTGCATTTGGGTTTGCATTTGGGTTT
// 1: CAAACCCAAATGCAAACCCAAATGCAAACCC:1 GGGTTTGCATTTGGGTTTGCATTTGGGTTTG
// ...
void db_nodes_print_verbose(const dBNode *nodes, size_t num,
                            const dBGraph *db_graph, FILE *out);

// Print in/outdegree - For debugging mostly
// indegree/outdegree (2 means >=2)
// 00: ! 01: + 02: {
// 10: - 11: = 12: <
// 20: } 21: > 22: *
void db_nodes_print_edges(const dBNode *nodes, size_t num,
                          const dBGraph *db_graph, FILE *out);


//
// Integrity checks
//

// Check an array of nodes denote a contigous path
// Prints warning and returns false on failure
bool db_node_check_nodes(const dBNode *nodes, size_t num,
                         const dBGraph *db_graph);

#endif /* DB_NODE_H_ */
