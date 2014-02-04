#include "global.h"
#include "binary_kmer.h"

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
boolean binary_kmer_less_than(BinaryKmer left, BinaryKmer right) {
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

// byte reverse complement look up table
// Example: byte representing ACTG -> CAGT
//   since A=00, C=01, G=10, T=11
//         ACTG = 00011110 (30)
//   revcmp[30] = 01001011 (0x4B)
//              => CAGT
static const uint8_t revcmp_table[256] =
{
  0xFF, 0xBF, 0x7F, 0x3F, 0xEF, 0xAF, 0x6F, 0x2F,
  0xDF, 0x9F, 0x5F, 0x1F, 0xCF, 0x8F, 0x4F, 0x0F,
  0xFB, 0xBB, 0x7B, 0x3B, 0xEB, 0xAB, 0x6B, 0x2B,
  0xDB, 0x9B, 0x5B, 0x1B, 0xCB, 0x8B, 0x4B, 0x0B,
  0xF7, 0xB7, 0x77, 0x37, 0xE7, 0xA7, 0x67, 0x27,
  0xD7, 0x97, 0x57, 0x17, 0xC7, 0x87, 0x47, 0x07,
  0xF3, 0xB3, 0x73, 0x33, 0xE3, 0xA3, 0x63, 0x23,
  0xD3, 0x93, 0x53, 0x13, 0xC3, 0x83, 0x43, 0x03,
  0xFE, 0xBE, 0x7E, 0x3E, 0xEE, 0xAE, 0x6E, 0x2E,
  0xDE, 0x9E, 0x5E, 0x1E, 0xCE, 0x8E, 0x4E, 0x0E,
  0xFA, 0xBA, 0x7A, 0x3A, 0xEA, 0xAA, 0x6A, 0x2A,
  0xDA, 0x9A, 0x5A, 0x1A, 0xCA, 0x8A, 0x4A, 0x0A,
  0xF6, 0xB6, 0x76, 0x36, 0xE6, 0xA6, 0x66, 0x26,
  0xD6, 0x96, 0x56, 0x16, 0xC6, 0x86, 0x46, 0x06,
  0xF2, 0xB2, 0x72, 0x32, 0xE2, 0xA2, 0x62, 0x22,
  0xD2, 0x92, 0x52, 0x12, 0xC2, 0x82, 0x42, 0x02,
  0xFD, 0xBD, 0x7D, 0x3D, 0xED, 0xAD, 0x6D, 0x2D,
  0xDD, 0x9D, 0x5D, 0x1D, 0xCD, 0x8D, 0x4D, 0x0D,
  0xF9, 0xB9, 0x79, 0x39, 0xE9, 0xA9, 0x69, 0x29,
  0xD9, 0x99, 0x59, 0x19, 0xC9, 0x89, 0x49, 0x09,
  0xF5, 0xB5, 0x75, 0x35, 0xE5, 0xA5, 0x65, 0x25,
  0xD5, 0x95, 0x55, 0x15, 0xC5, 0x85, 0x45, 0x05,
  0xF1, 0xB1, 0x71, 0x31, 0xE1, 0xA1, 0x61, 0x21,
  0xD1, 0x91, 0x51, 0x11, 0xC1, 0x81, 0x41, 0x01,
  0xFC, 0xBC, 0x7C, 0x3C, 0xEC, 0xAC, 0x6C, 0x2C,
  0xDC, 0x9C, 0x5C, 0x1C, 0xCC, 0x8C, 0x4C, 0x0C,
  0xF8, 0xB8, 0x78, 0x38, 0xE8, 0xA8, 0x68, 0x28,
  0xD8, 0x98, 0x58, 0x18, 0xC8, 0x88, 0x48, 0x08,
  0xF4, 0xB4, 0x74, 0x34, 0xE4, 0xA4, 0x64, 0x24,
  0xD4, 0x94, 0x54, 0x14, 0xC4, 0x84, 0x44, 0x04,
  0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20,
  0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00
};


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
  assert(seq != NULL);
  assert(strlen(seq) >= kmer_size);

  // Faster attempt
  size_t i;
  const char *k = seq, *end = seq + BKMER_TOP_BASES(kmer_size);
  BinaryKmer bkmer = BINARY_KMER_ZERO_MACRO;
  Nucleotide nuc;

  // Do first word
  for(; k < end; k++) {
    assert(char_is_acgt(*k));
    nuc = dna_char_to_nuc(*k);
    bkmer.b[0] = (bkmer.b[0] << 2) | nuc;
  }

  // Do remaining words
  for(i = 1; i < NUM_BKMER_WORDS; i++) {
    for(end += 32; k < end; k++) {
      assert(char_is_acgt(*k));
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
