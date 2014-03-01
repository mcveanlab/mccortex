#include "global.h"
#include "binary_kmer.h"
#include "binary_seq.h"

// bswap only used for approach 2 in binary_kmer_reverse_complement
#ifdef _MSC_VER
  #define bswap_32(x) _byteswap_ulong(x)
  #define bswap_64(x) _byteswap_uint64(x)
#elif defined(__APPLE__)
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

BinaryKmer binary_kmer_left_shift_add(const BinaryKmer bkmer, size_t kmer_size,
                                      Nucleotide nuc)
{
  BinaryKmer b = binary_kmer_left_shift_one_base(bkmer, kmer_size);
  b.b[NUM_BKMER_WORDS - 1] |= nuc;
  return b;
}

BinaryKmer binary_kmer_right_shift_add(const BinaryKmer bkmer, size_t kmer_size,
                                       Nucleotide nuc)
{
  BinaryKmer b = binary_kmer_right_shift_one_base(bkmer);
  b.b[0] |= ((uint64_t)(nuc)) << BKMER_TOP_BP_BYTEOFFSET(kmer_size);
  return b;
}

// uses revcmp_table in src/basic/binary_seq.c
BinaryKmer binary_kmer_reverse_complement(const BinaryKmer bkmer,
                                          size_t kmer_size)
{
  const size_t top_bits = BKMER_TOP_BITS(kmer_size), unused_bits = 64 - top_bits;
  size_t i, j, k;
  BinaryKmer revcmp = BINARY_KMER_ZERO_MACRO;
  uint64_t word;

  for(i = 0, j = NUM_BKMER_WORDS-1; i < NUM_BKMER_WORDS; i++, j--)
  {
    // swap word i into word j
    word = bkmer.b[i];

    // In word: reversing base order and bit negate (complementing base)
    // Approach 1 (using revcmp_table)
    for(k = 0; k < sizeof(uint64_t); k++)
    {
      revcmp.b[j] <<= 8;
      revcmp.b[j] |= revcmp_table[word & 0xff];
      word >>= 8;
    }
    // Approach 2
    // // Swap byte order
    // word = bswap_64(word);
    // // 4 bases within a byte, so swap their order
    // word = (((word & 0x0303030303030303UL) << 6) |
    //         ((word & 0x0c0c0c0c0c0c0c0cUL) << 2) |
    //         ((word & 0x3030303030303030UL) >> 2) |
    //         ((word & 0xc0c0c0c0c0c0c0c0UL) >> 6));
    // // Bitwise negate to complement bases
    // revcmp.b[j] = ~word;
  }

  // Now shift bits right by unused_bits
  for(i = NUM_BKMER_WORDS-1; i > 0; i--) {
    revcmp.b[i] = (revcmp.b[i] >> unused_bits) | (revcmp.b[i-1] << top_bits);
  }
  revcmp.b[0] >>= unused_bits;

  return revcmp;
}

// Get a random binary kmer -- useful for testing
BinaryKmer binary_kmer_random(size_t kmer_size)
{
  BinaryKmer bkmer = BINARY_KMER_ZERO_MACRO;
  size_t i;
  for(i = 0; i < sizeof(BinaryKmer) / 8; i++)
    bkmer.b[i] = (((uint64_t)rand()) << 32) | (uint64_t)rand();
  bkmer.b[0] >>= 64 - BKMER_TOP_BASES(kmer_size);
  return bkmer;
}

//
// Functions operating on strings
//

// Caller passes in preallocated BinaryKmer
// which is also returned in the return value
BinaryKmer binary_kmer_from_str(const char *seq, size_t kmer_size)
{
  ctx_assert(seq != NULL);
  ctx_assert(strlen(seq) >= kmer_size);

  // Faster attempt
  size_t i;
  const char *k = seq, *end = seq + BKMER_TOP_BASES(kmer_size);
  BinaryKmer bkmer = BINARY_KMER_ZERO_MACRO;
  Nucleotide nuc;

  // Do first word
  for(; k < end; k++) {
    ctx_assert(char_is_acgt(*k));
    nuc = dna_char_to_nuc(*k);
    bkmer.b[0] = (bkmer.b[0] << 2) | nuc;
  }

  // Do remaining words
  for(i = 1; i < NUM_BKMER_WORDS; i++) {
    for(end += 32; k < end; k++) {
      ctx_assert(char_is_acgt(*k));
      nuc = dna_char_to_nuc(*k);
      bkmer.b[i] = (bkmer.b[i] << 2) | nuc;
    }
  }

  return bkmer;
}

// Caller passes in allocated char* as 3rd argument which is then returned
// Note that the allocated space has to be kmer_size+1;
char *binary_kmer_to_str(const BinaryKmer bkmer, size_t kmer_size, char *seq)
{
  size_t i, j, k = kmer_size, topbases = BKMER_TOP_BASES(kmer_size);
  uint64_t word;

  // All but the top word
  for(i = NUM_BKMER_WORDS-1; i > 0; i--) {
    word = bkmer.b[i];
    for(j = 0; j < 32; j++) {
      seq[--k] = dna_nuc_to_char(word & 0x3);
      word >>= 2;
    }
  }

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
  size_t i, j, slen = (kmer_size+1)/2, k = slen;
  size_t toppairs = (BKMER_TOP_BASES(kmer_size)+1)/2;
  uint64_t word;

  for(i = NUM_BKMER_WORDS-1; i > 0; i--) {
    word = bkmer.b[i];
    for(j = 0; j < 32; j += 2) {
      seq[--k] = hex[word & 15];
      word >>= 4;
    }
  }

  // Top word
  word = bkmer.b[0];
  for(j = 0; j < toppairs; j++) {
    seq[--k] = hex[word & 15];
    word >>= 4;
  }

  seq[slen] = '\0';
}
