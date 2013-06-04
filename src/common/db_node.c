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
void bkmer_with_orientation(const BinaryKmer bkmer, Orientation orient,
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
// Coverages
//


void db_node_add_coverage(dBGraph *graph, hkey_t hkey, Colour col, long update)
{
  safe_add_covgs(graph->covgs[hkey][col], update, graph->covgs[hkey][col]);
}

void db_node_increment_coverage(dBGraph *graph, hkey_t hkey, Colour col)
{
  safe_covgs_increment(graph->covgs[hkey][col]);
}

Covg db_node_sum_covg_of_colours(const dBGraph *graph, hkey_t hkey,
                                 Colour first_col, int num_cols)
{
  Covg sum_covg = 0, *covgs = graph->covgs[hkey];
  Colour col, end;

  for(col = first_col, end = first_col+num_cols; col < end; col++)
    safe_add_covgs(sum_covg, covgs[col], sum_covg);

  return sum_covg;
}

Covg db_node_sum_covg_of_all_colours(const dBGraph *graph, hkey_t hkey)
{
  Covg sum_covg = 0, *covgs = graph->covgs[hkey];
  Colour col;

  for(col = 0; col < NUM_OF_COLOURS; col++)
    safe_add_covgs(sum_covg, covgs[col], sum_covg);

  return sum_covg;
}

Covg db_node_sum_covg_of_colourlist(const dBGraph *graph, hkey_t hkey,
                                    const Colour *cols, int num_cols)
{
  Covg sum_covg = 0, *covgs = graph->covgs[hkey];
  int i;

  for(i = 0; i < num_cols; i++)
    safe_add_covgs(sum_covg, covgs[cols[i]], sum_covg);

  return sum_covg;
}

//
// dBNodeBuffer
//
void db_node_buf_alloc(dBNodeBuffer*buf, size_t capacity)
{
  buf->capacity = ROUNDUP2POW(capacity);
  buf->nodes = malloc(sizeof(hkey_t) * buf->capacity);
  buf->orients = malloc(sizeof(Orientation) * buf->capacity);
  buf->len = 0;
  if(buf->nodes == NULL || buf->orients == NULL) die("Out of memory");
}

void db_node_buf_dealloc(dBNodeBuffer *buf)
{
  free(buf->nodes);
  free(buf->orients);
}

void db_node_buf_ensure_capacity(dBNodeBuffer *buf, size_t capacity)
{
  if(capacity > buf->capacity)
  {
    buf->capacity = ROUNDUP2POW(capacity);
    buf->nodes = realloc(buf->nodes, sizeof(hkey_t) * buf->capacity);
    buf->orients = realloc(buf->orients, sizeof(Orientation) * buf->capacity);
    if(buf->nodes == NULL || buf->orients == NULL) die("Out of memory");
  }
}

void db_node_buf_safe_add(dBNodeBuffer *buf, hkey_t node, Orientation orient)
{
  size_t n = buf->len;
  db_node_buf_ensure_capacity(buf, n+1);
  buf->nodes[n] = node;
  buf->orients[n] = orient;
  buf->len++;
}
