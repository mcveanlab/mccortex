#include "global.h"
#include "binary_kmer.h"

// This is exported
const BinaryKmer zero_bkmer = {0};

const char bnuc_to_char_array[4] = {'A','C','G','T'};

// 0:A, 1:C, 2:G, 3:T, 4:Undefined
const Nucleotide char_to_bnuc[128] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                      4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                      4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
                                      4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
                                      4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
                                      4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4};

#ifndef NDEBUG
// These have asserts in them
char binary_nuc_to_char(Nucleotide n)
{
  // Check only values 0..3 are being used
  assert((n & 0x3) == n);
  return bnuc_to_char_array[n];
}
#endif

// less than for 1 or 2 bitfields is defined in the header
#if NUM_BITFIELDS_IN_BKMER > 2
boolean binary_kmer_less_than(const BinaryKmer left, const BinaryKmer right)
{
  int i;

  // start at most significant end
  for(i = 0; i < NUM_BITFIELDS_IN_BKMER; i++)
  {
    if(left[i] < right[i]) return true;
    if(left[i] > right[i]) return false;
  }

  // if equal, return false ('less than' not 'less than or equal')
  return false;
}
#endif


void binary_kmer_right_shift_one_base(BinaryKmer kmer)
{
  int i;
  for(i = NUM_BITFIELDS_IN_BKMER - 1; i > 0; i--)
  {
    kmer[i] >>= 2;
    kmer[i] |= (kmer[i - 1] << 62);
  }

  kmer[0] >>= 2;
}

void binary_kmer_left_shift_one_base(BinaryKmer kmer, uint32_t kmer_size)
{
  int i;
  for(i = 0; i < NUM_BITFIELDS_IN_BKMER - 1; i++)
  {
    kmer[i] <<= 2;
    kmer[i] |= (kmer[i + 1] >> 62);
  }

  kmer[NUM_BITFIELDS_IN_BKMER - 1] <<= 2;

  // Mask top word
  kmer[0] &= (~(uint64_t)0 >> (64 - BKMER_TOP_BITS(kmer_size)));
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


// kmer and prealloc_revcmp_kmer may NOT point to the same address
BinaryKmerPtr binary_kmer_reverse_complement(const uint64_t *const restrict kmer,
                                             uint32_t kmer_size,
                                             BinaryKmerPtr restrict revcmp_kmer)
{
  assert(kmer != revcmp_kmer);
  // BinaryKmer kmer_copy;
  size_t i, j, k;

  for(i = 0; i < NUM_BITFIELDS_IN_BKMER; i++)
  {
    // swap word i into word j
    j = NUM_BITFIELDS_IN_BKMER - i - 1;
    uint64_t kmer_word = kmer[i];

    // Loop over bytes in binary kmer
    for(k = 0; k < sizeof(uint64_t); k++)
    {
      revcmp_kmer[j] <<= 8;
      revcmp_kmer[j] |= revcmp_table[kmer_word & 0xff];
      kmer_word >>= 8;
    }
  }

  // bits in top word
  size_t top_bits = BKMER_TOP_BITS(kmer_size);
  size_t unused_bits = 64 - top_bits;

  // Now shift bits right by unused_bits
  int w = NUM_BITFIELDS_IN_BKMER-1;

  revcmp_kmer[w] >>= unused_bits;

  for(w = w - 1; w >= 0; w--)
  {
    revcmp_kmer[w+1] |= revcmp_kmer[w] << top_bits;
    revcmp_kmer[w] >>= unused_bits;
  }

  return revcmp_kmer;
}

//
// Functions operating on strings
//

// Caller passes in preallocated BinaryKmer
// which is also returned in the return value
BinaryKmerPtr binary_kmer_from_str(const char *seq, uint32_t kmer_size,
                                   BinaryKmer prealloced_kmer)
{
  assert(seq != NULL);
  assert(prealloced_kmer != NULL);

  if(strlen(seq) < kmer_size) {
    printf("[kmer:%u>%zu]: '%s' %i\n", kmer_size, strlen(seq), seq, seq[0]);
  }

  assert(strlen(seq) >= kmer_size);

  uint32_t i;
  binary_kmer_init(prealloced_kmer);

  for(i = 0; i < kmer_size; i++)
  {
    Nucleotide nuc = binary_nuc_from_char(seq[i]);
    assert(nuc != Undefined);

    binary_kmer_left_shift_add(prealloced_kmer, kmer_size, nuc);
  }

  return prealloced_kmer;
}

// Caller passes in allocated char* as 3rd argument which is then returned
// User of this method is responsible for deallocating the returned sequence
// Note that the allocated space has to be kmer_size+1;
char *binary_kmer_to_str(const BinaryKmer bkmer, uint32_t kmer_size, char *seq)
{
  assert(seq != NULL);

  BinaryKmer local_bkmer;
  binary_kmer_assign(local_bkmer, bkmer);

  uint32_t i;
  for(i = kmer_size - 1; ; i--)
  {
    Nucleotide nuc = binary_kmer_last_nuc(local_bkmer);
    seq[i] = binary_nuc_to_char(nuc);
    binary_kmer_right_shift_one_base(local_bkmer);
    if(i == 0) break;
  }

  seq[kmer_size] = '\0';

  return seq;
}

void binary_nuc_from_str(Nucleotide *bases, const char *str, size_t len)
{
  size_t i;
  for(i = 0; i < len; i++)
    bases[i] = binary_nuc_from_char(str[i]);
}

void binary_nuc_to_str(const Nucleotide *bases, char *str, size_t len)
{
  size_t i;
  for(i = 0; i < len; i++)
    str[i] = binary_nuc_to_char(bases[i]);
}
