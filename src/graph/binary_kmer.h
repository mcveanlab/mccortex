#ifndef BINARY_KMER_H_
#define BINARY_KMER_H_

#include "bit_array/bit_macros.h"
#include "dna.h"

// coding is: [1]=xx665544 [0]=33221100

// NUM_BKMER_WORDS is the number of 64 bit words we use to encode a kmer
#define NUM_BKMER_WORDS (((MAX_KMER_SIZE)*2+63)/64)

typedef struct {
  uint64_t b[NUM_BKMER_WORDS];
} BinaryKmer;

// #define BKMER_BYTES (NUM_BKMER_WORDS * sizeof(uint64_t))
#define BKMER_BYTES (NUM_BKMER_WORDS * 8)

// BinaryKmer that is all zeros
extern const BinaryKmer zero_bkmer;

#define BINARY_KMER_ZERO_MACRO {.b = {0}}

// Hash functions
#if defined(USE_CITY_HASH)
  // Use Google's CityHash
  #include "misc/city.h"
  #define binary_kmer_hash(bkmer,rehash) (CityHash32((char*)bkmer.b, BKMER_BYTES) ^ rehash)
#elif defined(USE_XXHASH)
  // Use xxHash
  #include "xxHash/xxhash.h"
  #define binary_kmer_hash(bkmer,rehash) XXH32(bkmer.b, BKMER_BYTES, rehash)
#else
  // Use Bob Jenkin's lookup3
  #include "kmer_hash.h"
  #define binary_kmer_hash(bkmer,rehash) bklk3_hashlittle(bkmer, rehash)
#endif


// Since kmer_size is always odd, top word always has <= 62 bits used
// Number of bases store in all but the top word
#define BKMER_LOWER_BASES              ((NUM_BKMER_WORDS-1)*32)
#define BKMER_TOP_BASES(ksize)         ((ksize)&31)
#define BKMER_TOP_BITS(ksize)          (BKMER_TOP_BASES(ksize) * 2)
#define BKMER_TOP_BP_BYTEOFFSET(ksize) (BKMER_TOP_BITS(ksize) - 2)

#define binary_kmer_first_nuc(bkmer,ksize) \
  (((bkmer).b[0] >> BKMER_TOP_BP_BYTEOFFSET(ksize)) & 0x3)

#define binary_kmer_last_nuc(bkmer)  ((bkmer).b[NUM_BKMER_WORDS - 1] & 0x3)

#define binary_kmer_set_first_nuc(bkmer,nuc,ksize)                              \
  ((bkmer)->b[0] = ((bkmer)->b[0] & bitmask64(BKMER_TOP_BP_BYTEOFFSET(ksize))) |\
                   (((uint64_t)(nuc)) << BKMER_TOP_BP_BYTEOFFSET(ksize)))

#define binary_kmer_set_last_nuc(bkmer,nuc) \
        ((bkmer)->b[NUM_BKMER_WORDS - 1] \
           = ((bkmer)->b[NUM_BKMER_WORDS - 1] & 0xfffffffffffffffcUL) | (nuc))

#if NUM_BKMER_WORDS == 1
  #define binary_kmers_cmp(x,y) ((x).b[0] < (y).b[0] ? -1 : (x).b[0] > (y).b[0])
#else
  int binary_kmers_cmp(BinaryKmer a, BinaryKmer b);
#endif

#if NUM_BKMER_WORDS == 1
  #define binary_kmers_are_equal(x,y) ((x).b[0] == (y).b[0])
  #define binary_kmer_is_zero(x)      ((x).b[0] == 0UL)
  #define binary_kmer_less_than(x,y)  ((x).b[0] < (y).b[0])
#elif NUM_BKMER_WORDS == 2
  #define binary_kmers_are_equal(x,y) ((x).b[0]==(y).b[0] && (x).b[1]==(y).b[1])
  #define binary_kmer_is_zero(x)      (((x).b[0] | (x).b[1]) == 0UL)
  #define binary_kmer_less_than(x,y) \
          ((x).b[0] < (y).b[0] || ((x).b[0] == (y).b[0] && (x).b[1] < (y).b[1]))
#else /* NUM_BKMER_WORDS > 2 */
  #define binary_kmers_are_equal(x,y) (memcmp((x).b,(y).b,BKMER_BYTES) == 0)
  #define binary_kmer_is_zero(x)      binary_kmers_are_equal((x), zero_bkmer)
  bool binary_kmer_less_than(BinaryKmer left, BinaryKmer right);
#endif

#define binary_kmer_oversized(bk,k)  ((bk).b[0] & (UINT64_MAX << BKMER_TOP_BITS(k)))

static inline int binary_kmers_qcmp(const void *aa, const void *bb)
{
  const BinaryKmer *a = (const BinaryKmer*)aa, *b = (const BinaryKmer*)bb;
  return binary_kmers_cmp(*a, *b);
}

//
// Functions
//

// For a given kmer, get the BinaryKmer 'key':
// the lower of the kmer vs reverse complement of itself
BinaryKmer binary_kmer_get_key(const BinaryKmer kmer, size_t kmer_size);

#if NUM_BKMER_WORDS == 1

  static inline BinaryKmer binary_kmer_right_shift_one_base(const BinaryKmer bkmer)
  {
    BinaryKmer b = bkmer; b.b[0] >>= 2; return b;
  }

  static inline BinaryKmer binary_kmer_left_shift_one_base(const BinaryKmer bkmer,
                                                           size_t kmer_size)
  {
    BinaryKmer b = bkmer;
    b.b[NUM_BKMER_WORDS - 1] <<= 2;
    b.b[0] &= (UINT64_MAX >> (64 - BKMER_TOP_BITS(kmer_size)));
    return b;
  }

#else

  // CTAGT -> ACTAG (add blank 'A' to first position)
  // Shift towards most significant position
  BinaryKmer binary_kmer_right_shift_one_base(const BinaryKmer bkmer);

  // CTAGT -> TAGTA (add blank 'A' to last position)
  BinaryKmer binary_kmer_left_shift_one_base(const BinaryKmer bkmer,
                                             size_t kmer_size);

#endif

static inline
BinaryKmer binary_kmer_left_shift_add(const BinaryKmer bkmer, size_t kmer_size,
                                      Nucleotide nuc)
{
  BinaryKmer b = binary_kmer_left_shift_one_base(bkmer, kmer_size);
  b.b[NUM_BKMER_WORDS - 1] |= nuc;
  return b;
}

static inline
BinaryKmer binary_kmer_right_shift_add(const BinaryKmer bkmer, size_t kmer_size,
                                       Nucleotide nuc)
{
  BinaryKmer b = binary_kmer_right_shift_one_base(bkmer);
  b.b[0] |= ((uint64_t)(nuc)) << BKMER_TOP_BP_BYTEOFFSET(kmer_size);
  return b;
}

// Reverse complement a binary kmer from kmer into revcmp_kmer
BinaryKmer binary_kmer_reverse_complement(const BinaryKmer bkmer, size_t kmer_size);

// Get a random binary kmer -- useful for testing
BinaryKmer binary_kmer_random(size_t kmer_size);

// BinaryKmer <-> String functions
char* binary_kmer_to_str(const BinaryKmer kmer, size_t kmer_size, char *seq);
BinaryKmer binary_kmer_from_str(const char *seq, size_t kmer_size);

void binary_kmer_to_hex(const BinaryKmer bkmer, size_t kmer_size, char *seq);

#endif /* BINARY_KMER_H_ */
