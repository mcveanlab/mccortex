#include "global.h"
#include "db_node.h"

// Only print covg overflow warning once
static char overflow_warning_printed = 0;

#define print_overflow_warning() do {                                          \
  if(!overflow_warning_printed) {                                              \
    warn("caught integer overflow (some coverages may be underestimates)");    \
    overflow_warning_printed = 1;                                              \
  }                                                                            \
} while(0)

#define safe_add_covgs(covg1,covg2,sum) do {                                   \
  if(COVG_MAX - (covg1) >= (covg2)) { (sum) = (covg1) + (covg2); }             \
  else { (sum) = COVG_MAX; print_overflow_warning(); }                         \
} while(0)

#define safe_covgs_increment(covg) do {                                        \
  if((covg) == COVG_MAX) { print_overflow_warning(); }                         \
  else (covg)++;                                                               \
} while(0)


// For a given kmer, get the BinaryKmer 'key':
// the lower of the kmer vs reverse complement of itself
// kmer and kmer_key must NOT point to overlapping memory
Key db_node_get_key(const uint64_t* const restrict kmer, uint32_t kmer_size,
                    Key restrict kmer_key)
{
  assert(kmer != kmer_key);

  /*
  // Simpler version is slightly slower
  binary_kmer_reverse_complement(kmer, kmer_size, kmer_key);
  if(binary_kmer_less_than(kmer, kmer_key))
    binary_kmer_assign(kmer_key, kmer);
  */

  // Get first and last nucleotides
  Nucleotide first = binary_kmer_first_nuc(kmer, kmer_size);
  Nucleotide last = binary_kmer_last_nuc(kmer);
  Nucleotide rev_last = binary_nuc_complement(last);

  if(first < rev_last) {
    binary_kmer_assign(kmer_key, kmer);
  }
  // else if(first > rev_last) { // reduce conditionals to speed up
  //   binary_kmer_reverse_complement(kmer, kmer_size, kmer_key);
  // }
  else
  {
    // Don't know which is going to be correct -- this will happen 1 in 4 times
    binary_kmer_reverse_complement(kmer, kmer_size, kmer_key);
    if(binary_kmer_less_than(kmer, kmer_key)) {
      binary_kmer_assign(kmer_key, kmer);
    }
  }

  return kmer_key;
}

//
// Orientations
//
void db_node_oriented_bkmer(const BinaryKmer bkmer, Orientation orient,
                            uint32_t kmer_size, BinaryKmer result)
{
  if(orient == forward) {
    binary_kmer_assign(result, bkmer);
  }
  else {
    binary_kmer_reverse_complement(bkmer, kmer_size, result);
  }
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
  *nucleotide = __builtin_ctz(edges);

  return (b == 0);
}

//
// Edges per colour
//

Edges edges_get_union(const Edges *edges_arr, size_t num)
{
  // Slow version
  // Edges edges = 0;
  // size_t i;
  // for(i = 0; i < NUM_OF_COLOURS; i++) edges |= edges_arr[i];
  // return edges;

  // Faster version deals with 8 bytes at a time
  uint64_t i, edges = 0;
  const uint64_t *ptr = (const uint64_t*)edges_arr;

  for(i = 0; i < num / 8; i++, ptr++) edges |= *ptr;

  edges |= *ptr & (~(uint64_t)0 >> (64-(8*(num % 8))));
  edges |= edges >> 32;
  edges |= edges >> 16;
  edges |= edges >> 8;
  return edges & 0xFF;
}

//
// Coverages
//


void db_node_add_coverage(dBGraph *graph, hkey_t hkey, Colour col, long update)
{
  safe_add_covgs(db_node_col_covg(graph,hkey,col), update,
                 db_node_col_covg(graph,hkey,col));
}

void db_node_increment_coverage(dBGraph *graph, hkey_t hkey, Colour col)
{
  safe_covgs_increment(db_node_col_covg(graph,hkey,col));
}

Covg db_node_sum_covg_of_colours(const dBGraph *graph, hkey_t hkey,
                                 Colour first_col, int num_of_cols)
{
  const Covg *covgs = &db_node_col_covg(graph,hkey,0);

  Covg sum_covg = 0;
  Colour col, end;

  for(col = first_col, end = first_col+num_of_cols; col < end; col++)
    safe_add_covgs(sum_covg, covgs[col], sum_covg);

  return sum_covg;
}


//
// dBNodeBuffer
//
void db_node_buf_alloc(dBNodeBuffer *buf, size_t capacity)
{
  buf->capacity = ROUNDUP2POW(capacity);
  buf->data = malloc(sizeof(dBNode) * buf->capacity);
  buf->len = 0;
  if(buf->data == NULL) die("Out of memory");
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
    buf->data = realloc(buf->data, sizeof(dBNode) * buf->capacity);
    if(buf->data == NULL) die("Out of memory");
  }
}

void db_node_buf_safe_add(dBNodeBuffer *buf, hkey_t node, Orientation orient)
{
  size_t n = buf->len;
  db_node_buf_ensure_capacity(buf, n+1);
  buf->data[n].node = node;
  buf->data[n].orient = orient;
  buf->len++;
}

void db_nodes_to_str(const dBNode *nodes, size_t num,
                     const dBGraph *db_graph, char *str)
{
  if(num == 0) return;

  size_t i;
  uint32_t kmer_size = db_graph->kmer_size;
  ConstBinaryKmerPtr bkmer = db_node_bkmer(db_graph, nodes[0].node);
  Nucleotide nuc;

  binary_kmer_to_str(bkmer, kmer_size, str);
  if(nodes[0].orient == reverse) reverse_complement_str(str, kmer_size);

  for(i = 1; i < num; i++) {
    bkmer = db_node_bkmer(db_graph, nodes[i].node);
    nuc = db_node_last_nuc(bkmer, nodes[i].orient, kmer_size);
    str[kmer_size+i-1] = binary_nuc_to_char(nuc);
  }

  str[kmer_size+num-1] = '\0';
}
