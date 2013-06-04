#ifndef BINARY_KMER_H_
#define BINARY_KMER_H_

#include <inttypes.h>

// Think of NUM_BITFIELDS_IN_BKMER as the number of uint64_t we
// encode the kmer in
typedef uint64_t BinaryKmer[NUM_BITFIELDS_IN_BKMER];
typedef uint64_t* BinaryKmerPtr;

// BinaryKmer that is all zeros
extern const BinaryKmer zero_bkmer;

// Encoding
typedef enum
{
  Adenine   = 0,
  Cytosine  = 1,
  Guanine   = 2,
  Thymine   = 3,
  Undefined = 4,
} Nucleotide;

// Hash functions
#ifdef CITY_HASH
  // Use Google's CityHash
  #include "city.h"
  #define binary_kmer_hash(key,rehash,bit_mask) \
    (CityHash32((char*)key, sizeof(BinaryKmer)) ^ (rehash)) & (bit_mask))
#else
  // Use Bob Jenkin's lookup3
  #include "lookup3.h"
  #define binary_kmer_hash(key,rehash,bit_mask) \
    (hashlittle(key, sizeof(BinaryKmer), (rehash)) & (bit_mask))
#endif

extern const char bnuc_to_char_array[4];
extern const Nucleotide char_to_bnuc[128];
extern const char complement_base[128];

#define is_base_char(c) (char_to_bnuc[(int)(c)] != Undefined)
#define binary_nucleotide_complement(n) (~(n) & 0x3)
#define char_to_binary_nucleotide(c)  char_to_bnuc[(int)(c)]

#ifdef NDEBUG
  #define binary_nucleotide_to_char(c)  bnuc_to_char_array[(c)]
  #define char_nucleotide_complement(c) complement_base[(int)(c)]
#else
  char binary_nucleotide_to_char(Nucleotide n);
  char char_nucleotide_complement(char c);
#endif

// Number of bases store in all but the top word
#define BKMER_LOWER_BASES              ((NUM_BITFIELDS_IN_BKMER-1)*32)
#define BKMER_TOP_BASES(ksize)         ((ksize) - BKMER_LOWER_BASES)
#define BKMER_TOP_BITS(ksize)          (BKMER_TOP_BASES(ksize) << 1)
#define BKMER_TOP_BP_BYTEOFFSET(ksize) (BKMER_TOP_BITS(ksize) - 2)

#define binary_kmer_init(b)          memset((b), 0, sizeof(BinaryKmer))
#define binary_kmer_assign(tgt,src)  memcpy((tgt), (src), sizeof(BinaryKmer))
#define binary_kmers_cmp(a,b)        memcmp((a),(b),sizeof(BinaryKmer))

#define binary_kmer_first_nucleotide(bkmer,ksize) \
        (((bkmer)[0] >> BKMER_TOP_BP_BYTEOFFSET(ksize)) & 0x3)
#define binary_kmer_last_nuc(bkmer) \
        ((bkmer)[NUM_BITFIELDS_IN_BKMER - 1] & 0x3)

#define binary_kmer_set_first_nucleotide(bkmer,nuc,ksize) \
        ((bkmer)[0] |= (nuc) << BKMER_TOP_BP_BYTEOFFSET(ksize))
#define binary_kmer_set_last_nuc(bkmer,nuc) \
        ((bkmer)[NUM_BITFIELDS_IN_BKMER - 1] |= (nuc))

#if NUM_BITFIELDS_IN_BKMER == 1
  #define binary_kmers_are_equal(a,b) ((a)[0] == (b)[0])
  #define binary_kmer_is_zero(b)      ((b)[0] == 0UL)
  #define binary_kmer_less_than(a,b)  ((a)[0] < (b)[0])
#elif NUM_BITFIELDS_IN_BKMER == 2
  #define binary_kmers_are_equal(a,b) ((a)[0] == (b)[0] && (a)[1] == (b)[1])
  #define binary_kmer_is_zero(b)      (((b)[0] | (b)[1]) == 0UL)
  #define binary_kmer_less_than(a,b) \
          ((a)[0] < (b)[0] || ((a)[0] == (b)[0] && (a)[1] < (b)[1]))
#else /* NUM_BITFIELDS_IN_BKMER > 2 */
  #define binary_kmers_are_equal(a,b) (binary_kmers_cmp((a),(b)) == 0)
  #define binary_kmer_is_zero(b)      binary_kmers_are_equal((b), zero_bkmer)
  boolean binary_kmer_less_than(const BinaryKmer left, const BinaryKmer right);
#endif

// Functions
void binary_kmer_right_shift_one_base(BinaryKmer kmer);
void binary_kmer_left_shift_one_base(BinaryKmer kmer, uint32_t kmer_size);

// Reverse complement a binary kmer from kmer into revcmp_kmer
// kmer and revcmp_kmer must NOT point to the same address
BinaryKmerPtr binary_kmer_reverse_complement(const uint64_t *const restrict kmer,
                                             uint32_t kmer_size,
                                             BinaryKmerPtr restrict revcmp_kmer);

// BinaryKmer <-> String functions
char *binary_kmer_to_str(const BinaryKmer kmer, uint32_t kmer_size, char *seq);
BinaryKmerPtr binary_kmer_from_str(const char *seq, uint32_t kmer_size,
                                   BinaryKmer prealloced_kmer);

char *reverse_complement_str(char *str, size_t length);

#endif /* BINARY_KMER_H_ */
