#include "global.h"
#include "binary_kmer.h"
#include "binary_seq.h"

#if defined(__APPLE__)
  #include <libkern/OSByteOrder.h>
  #define bswap_32(x) OSSwapInt32(x)
  #define bswap_64(x) OSSwapInt64(x)
#else
  #include <byteswap.h>
#endif

// This is exported
const BinaryKmer zero_bkmer = BINARY_KMER_ZERO_MACRO;

// less than for 1 or 2 bitfields is defined in the header
#if NUM_BKMER_WORDS > 2
bool binary_kmer_less_than(BinaryKmer left, BinaryKmer right) {
  size_t i;
  for(i = 0; i < NUM_BKMER_WORDS && left.b[i] == right.b[i]; i++);
  return (i < NUM_BKMER_WORDS && left.b[i] < right.b[i]);
}
#endif


BinaryKmer binary_kmer_from_old(BinaryKmer bkmer, size_t kmer_size)
{
  size_t o = 0, x = 2*(kmer_size&31);
  BinaryKmer nbkmer = {{0}};
  for(o = 0; o < x; o+=2) nbkmer.b[0] |= ((bkmer.b[0]>>(x-o-2))&3) << o;
  #if NUM_BKMER_WORDS > 1
  size_t i, j, w = 0;
  for(i = 0; i < NUM_BKMER_WORDS-1; i++) {
    for(j = 0; j < 64; j+=2) {
      nbkmer.b[w] |= ((bkmer.b[i]>>j)&3) << o;
      o += 2; w += (o & 64); o = (o & 63);
    }
  }
  #endif
  return nbkmer;
}

BinaryKmer binary_kmer_to_old(BinaryKmer bkmer, size_t kmer_size)
{
  BinaryKmer nbkmer = {{0}};
  (void)bkmer; (void)kmer_size;
  return nbkmer;
}

#if NUM_BKMER_WORDS > 1
  int binary_kmers_cmp(BinaryKmer a, BinaryKmer b)
  {
    size_t i;
    for(i = 0; i < NUM_BKMER_WORDS; i++) {
      if(a.b[i] < b.b[i]) return -1;
      if(a.b[i] > b.b[i]) return  1;
    }
    return 0;
  }
#endif

// For a given kmer, get the BinaryKmer 'key':
// the lower of the kmer vs reverse complement of itself
// kmer and kmer_key must NOT point to overlapping memory
BinaryKmer binary_kmer_get_key(const BinaryKmer bkmer, size_t kmer_size)
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

#if NUM_BKMER_WORDS > 1

// CTAGT -> ACTAG (add blank 'A' to first position)
// Shift towards most significant position
BinaryKmer binary_kmer_right_shift_one_base(const BinaryKmer bkmer)
{
  BinaryKmer b = bkmer;
  size_t i;

  for(i = NUM_BKMER_WORDS - 1; i > 0; i--)
  {
    b.b[i] >>= 2;
    b.b[i] |= (b.b[i - 1] << 62);
  }

  b.b[0] >>= 2;
  return b;
}

// CTAGT -> TAGTA (add blank 'A' to last position)
// Shift towards least significant position
BinaryKmer binary_kmer_left_shift_one_base(const BinaryKmer bkmer,
                                           size_t kmer_size)
{
  BinaryKmer b = bkmer;
  size_t i;

  for(i = 0; i+1 < NUM_BKMER_WORDS; i++)
  {
    b.b[i] <<= 2;
    b.b[i] |= (b.b[i + 1] >> 62);
  }

  b.b[NUM_BKMER_WORDS - 1] <<= 2;

  // Mask top word
  b.b[0] &= (~(uint64_t)0 >> (64 - BKMER_TOP_BITS(kmer_size)));
  return b;
}

#endif /* NUM_BKMER_WORDS > 1 */

// For profiling see dev/bkmer_revcmp/
BinaryKmer binary_kmer_reverse_complement(const BinaryKmer bkmer,
                                          size_t kmer_size)
{
  const size_t top_bits = BKMER_TOP_BITS(kmer_size), unused_bits = 64 - top_bits;
  size_t i, j;
  BinaryKmer revcmp = BINARY_KMER_ZERO_MACRO;
  uint64_t word;

  for(i = 0, j = NUM_BKMER_WORDS-1; i < NUM_BKMER_WORDS; i++, j--)
  {
    // Swap byte order
    word = bswap_64(bkmer.b[i]);
    // 4 bases within a byte, so swap their order
    word = (((word & 0x0303030303030303UL) << 6) |
            ((word & 0x0c0c0c0c0c0c0c0cUL) << 2) |
            ((word & 0x3030303030303030UL) >> 2) |
            ((word & 0xc0c0c0c0c0c0c0c0UL) >> 6));
    // Bitwise negate to complement bases
    revcmp.b[j] = ~word;
  }

#if NUM_BKMER_WORDS > 1
  // Need to shift right
  for(i = NUM_BKMER_WORDS-1; i > 0; i--) {
    revcmp.b[i] = (revcmp.b[i] >> unused_bits) | (revcmp.b[i-1] << top_bits);
  }
#endif

  revcmp.b[0] >>= unused_bits;

  return revcmp;
}

// Get a random binary kmer -- useful for testing
BinaryKmer binary_kmer_random(size_t kmer_size)
{
  BinaryKmer bkmer = BINARY_KMER_ZERO_MACRO;
  size_t i;
  for(i = 0; i < NUM_BKMER_WORDS; i++)
    bkmer.b[i] = (((uint64_t)rand()) << 32) | (uint64_t)rand();
  bkmer.b[0] >>= 64 - BKMER_TOP_BASES(kmer_size);
  return bkmer;
}

//
// Functions operating on strings
//

BinaryKmer binary_kmer_from_str(const char *seq, size_t kmer_size)
{
  ctx_assert(seq != NULL);
  ctx_assert(strlen(seq) >= kmer_size);

  // Faster attempt
  const char *k = seq, *end = seq + BKMER_TOP_BASES(kmer_size);
  BinaryKmer bkmer = BINARY_KMER_ZERO_MACRO;
  Nucleotide nuc;

  // Do first (partial) word
  for(; k < end; k++) {
    ctx_assert(char_is_acgt(*k));
    nuc = dna_char_to_nuc(*k);
    bkmer.b[0] = (bkmer.b[0] << 2) | nuc;
  }

#if NUM_BKMER_WORDS > 1
  // Do remaining (whole) words
  size_t i;
  for(i = 1; i < NUM_BKMER_WORDS; i++) {
    for(end += 32; k < end; k++) {
      ctx_assert(char_is_acgt(*k));
      nuc = dna_char_to_nuc(*k);
      bkmer.b[i] = (bkmer.b[i] << 2) | nuc;
    }
  }
#endif

  return bkmer;
}

// Caller passes in allocated char* as 3rd argument which is then returned
// Note that the allocated space has to be kmer_size+1;
char *binary_kmer_to_str(const BinaryKmer bkmer, size_t kmer_size, char *seq)
{
  size_t j, k = kmer_size, topbases = BKMER_TOP_BASES(kmer_size);
  uint64_t word;

#if NUM_BKMER_WORDS > 1
  // All but the top word
  size_t i;
  for(i = NUM_BKMER_WORDS-1; i > 0; i--) {
    word = bkmer.b[i];
    for(j = 0; j < 32; j++) {
      seq[--k] = dna_nuc_to_char(word & 0x3);
      word >>= 2;
    }
  }
#endif

  // Top word
  word = bkmer.b[0];
  for(j = 0; j < topbases; j++) {
    seq[--k] = dna_nuc_to_char(word & 0x3);
    word >>= 2;
  }

  seq[kmer_size] = '\0';

  return seq;
}

static const char hex[16] = "0123456789abcdef";

void binary_kmer_to_hex(const BinaryKmer bkmer, size_t kmer_size, char *seq)
{
  size_t j, slen = (kmer_size+1)/2, k = slen;
  size_t toppairs = (BKMER_TOP_BASES(kmer_size)+1)/2;
  uint64_t word;

#if NUM_BKMER_WORDS > 1
  size_t i;
  for(i = NUM_BKMER_WORDS-1; i > 0; i--) {
    word = bkmer.b[i];
    for(j = 0; j < 32; j += 2) {
      seq[--k] = hex[word & 15];
      word >>= 4;
    }
  }
#endif

  // Top word
  word = bkmer.b[0];
  for(j = 0; j < toppairs; j++) {
    seq[--k] = hex[word & 15];
    word >>= 4;
  }

  seq[slen] = '\0';
}
