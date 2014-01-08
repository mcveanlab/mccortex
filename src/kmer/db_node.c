#include "global.h"
#include <ctype.h> // toupper
#include "db_graph.h"
#include "db_node.h"
#include "util.h"

// For a given kmer, get the BinaryKmer 'key':
// the lower of the kmer vs reverse complement of itself
// kmer and kmer_key must NOT point to overlapping memory
BinaryKmer db_node_get_key(const BinaryKmer bkmer, size_t kmer_size)
{
  BinaryKmer bkey;

  // Get first and last nucleotides
  Nucleotide first = binary_kmer_first_nuc(bkmer, kmer_size);
  Nucleotide last = binary_kmer_last_nuc(bkmer);
  Nucleotide rev_last = dna_nuc_complement(last);

  if(first < rev_last) return bkmer;

  // Don't know which is going to be correct -- this will happen 1 in 4 times
  bkey = binary_kmer_reverse_complement(bkmer, kmer_size);
  return (binary_kmer_less_than(bkmer, bkey) ? bkmer : bkey);
}

//
// Edges
//

boolean edges_has_precisely_one_edge(Edges edges, Orientation orientation,
                                     Nucleotide *nucleotide)
{
  edges = edges_with_orientation(edges, orientation);
  if(edges == 0) return false;

  // bit-hack: b == 0 iff number of edges is 1
  // 1100 & 1011 => 1000
  // 0100 & 0011 => 0000
  int b = edges & (edges - 1);

  // __builtin_ctz returns number of trailing zeros
  // 0b0001 returns 0
  // 0b0010 returns 1
  // 0b0100 returns 2
  // 0b1000 returns 3
  *nucleotide = (Nucleotide)__builtin_ctz(edges);

  return (b == 0);
}

//
// Edges per colour
//

Edges edges_get_union(const Edges *edges_arr, size_t num)
{
  if(num == 1) return edges_arr[0];

  // Slow version
  // Edges edges = 0; size_t i;
  // for(i = 0; i < NUM_OF_COLOURS; i++) edges |= edges_arr[i];
  // return edges;

  // Faster version deals with 8 bytes at a time
  uint64_t i, edges = 0, tmp, end = 8*(num/8);

  for(i = 0; i < end; i += 8) {
    memcpy(&tmp, edges_arr+i, sizeof(uint64_t));
    edges |= tmp;
  }

  memcpy(&tmp, edges_arr+i, sizeof(uint64_t));
  edges |= tmp & (~(uint64_t)0 >> (64-(8*(num % 8))));

  // with unaligned memory access
  // const uint64_t *ptr = (const uint64_t*)((size_t)edges_arr);
  // for(i = 0; i < num / 8; i++, ptr++) edges |= *ptr;
  // edges |= *ptr & (~(uint64_t)0 >> (64-(8*(num % 8))));

  edges |= edges >> 32;
  edges |= edges >> 16;
  edges |= edges >> 8;
  return edges & 0xFF;
}

// kmer_col_edge_str should be 9 chars long
// Return pointer to kmer_col_edge_str
char* db_node_get_edges_str(Edges edges, char* kmer_col_edge_str)
{
  int i;
  char str[] = "acgt";

  char left = edges >> 4;
  left = (char)rev_nibble(left);
  char right = edges & 0xf;

  for(i = 0; i < 4; i++)
    kmer_col_edge_str[i] = (left & (0x1 << i) ? str[i] : '.');

  for(i = 0; i < 4; i++)
    kmer_col_edge_str[i+4] = (char)toupper(right & (0x1 << i) ? str[i] : '.');

  kmer_col_edge_str[8] = '\0';

  return kmer_col_edge_str;
}

//
// Coverages
//

static inline void safe_add_covg(Covg *a, Covg b)
{
  uint64_t y = (uint64_t)*a+b;
  *a = y > COVG_MAX ? COVG_MAX : (Covg)y;
}

void db_node_add_col_covg(dBGraph *graph, hkey_t hkey, Colour col, Covg update)
{
  safe_add_covg(&db_node_col_covg(graph,col,hkey), update);
}

void db_node_increment_coverage(dBGraph *graph, hkey_t hkey, Colour col)
{
  safe_add_covg(&db_node_col_covg(graph,col,hkey), 1);
}

Covg db_node_sum_covg(const dBGraph *graph, hkey_t hkey)
{
  const Covg *covgs = &db_node_col_covg(graph,0,hkey);
  Covg sum_covg = 0;
  size_t col;

  for(col = 0; col < graph->num_of_cols; col++)
    safe_add_covg(&sum_covg, covgs[col]);

  return sum_covg;
}

//
// dBNodeBuffer
//
void db_node_buf_alloc(dBNodeBuffer *buf, size_t capacity)
{
  buf->capacity = ROUNDUP2POW(capacity);
  buf->data = malloc2(sizeof(dBNode) * buf->capacity);
  buf->len = 0;
}

void db_node_buf_dealloc(dBNodeBuffer *buf)
{
  free(buf->data);
}

void db_node_buf_ensure_capacity(dBNodeBuffer *buf, size_t capacity)
{
  if(capacity > buf->capacity)
  {
    buf->capacity = ROUNDUP2POW(capacity);
    buf->data = realloc2(buf->data, sizeof(dBNode) * buf->capacity);
  }
}

void db_node_buf_safe_add(dBNodeBuffer *buf, hkey_t node, Orientation orient)
{
  size_t n = buf->len;
  db_node_buf_ensure_capacity(buf, n+1);
  buf->data[n].key = node;
  buf->data[n].orient = orient;
  buf->len++;
}

void db_nodes_to_str(const dBNode *nodes, size_t num,
                     const dBGraph *db_graph, char *str)
{
  if(num == 0) return;

  size_t i;
  size_t kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer = db_node_bkmer(db_graph, nodes[0].key);
  Nucleotide nuc;

  binary_kmer_to_str(bkmer, kmer_size, str);
  if(nodes[0].orient == REVERSE) dna_reverse_complement_str(str, kmer_size);

  for(i = 1; i < num; i++) {
    bkmer = db_node_bkmer(db_graph, nodes[i].key);
    nuc = db_node_last_nuc(bkmer, nodes[i].orient, kmer_size);
    str[kmer_size+i-1] = dna_nuc_to_char(nuc);
  }

  str[kmer_size+num-1] = '\0';
}

void db_nodes_print(const dBNode *nodes, size_t num,
                    const dBGraph *db_graph, FILE *out)
{
  const size_t kmer_size = db_graph->kmer_size;
  size_t i;
  Nucleotide nuc;
  BinaryKmer bkmer;
  char tmp[MAX_KMER_SIZE+1];

  bkmer = db_graph_oriented_bkmer(db_graph, nodes[0].key, nodes[0].orient);
  binary_kmer_to_str(bkmer, kmer_size, tmp);
  fputs(tmp, out);

  for(i = 1; i < num; i++) {
    bkmer = db_node_bkmer(db_graph, nodes[i].key);
    nuc = db_node_last_nuc(bkmer, nodes[i].orient, kmer_size);
    fputc(dna_nuc_to_char(nuc), out);
  }
}

void db_nodes_gzprint(const dBNode *nodes, size_t num,
                      const dBGraph *db_graph, gzFile out)
{
  size_t i, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkmer;
  char tmp[MAX_KMER_SIZE+1];

  bkmer = db_graph_oriented_bkmer(db_graph, nodes[0].key, nodes[0].orient);
  binary_kmer_to_str(bkmer, kmer_size, tmp);
  gzputs(out, tmp);

  for(i = 1; i < num; i++) {
    bkmer = db_node_bkmer(db_graph, nodes[i].key);
    nuc = db_node_last_nuc(bkmer, nodes[i].orient, kmer_size);
    gzputc(out, dna_nuc_to_char(nuc));
  }
}
