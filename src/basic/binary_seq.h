#ifndef BINARY_SEQ_H_
#define BINARY_SEQ_H_

#include "dna.h"

extern const uint8_t revcmp_table[256];

// assume nbases > 0
#define bases_in_top_byte(nbases) ((((nbases) - 1) & 3) + 1)
#define bits_in_top_byte(nbases) (bases_in_top_byte(nbases) * 2)

#define binary_seq_unpack_byte(arr,b) do {   \
  (arr)[0] =  (b)     & 0x3;                 \
  (arr)[1] = ((b)>>2) & 0x3;                 \
  (arr)[2] = ((b)>>4) & 0x3;                 \
  (arr)[3] =  (b)>>6;                        \
} while(0)

#define binary_seq_pack_byte(arr) \
        (((arr)[3]<<6) | ((arr)[2]<<4) | ((arr)[1]<<2) | ((arr)[0]))

// Fetch a given base. Four bases per byte
// ptr should point directly to sequence
static inline Nucleotide binary_seq_get(const uint8_t *ptr, size_t idx)
{
  size_t byte = idx / 4, offset = (idx & 3)*2;
  return (ptr[byte] >> offset) & 3;
}

// Add a given base. Four bases per byte
// ptr should point directly to sequence
// Doesn't provide any masking - should already be zeroed
static inline void binary_seq_set(uint8_t *ptr, size_t idx, Nucleotide nuc)
{
  ctx_assert(nuc < 4);
  size_t byte = idx / 4, offset = (idx & 3)*2;
  // 11111100 11110011 11001111 00111111
  const uint8_t masks[4] = {0xfc, 0xf3, 0xcf, 0x3f};
  ptr[byte] = (ptr[byte] & masks[idx&3]) | (nuc << offset);
}

void binary_seq_reverse_complement(uint8_t *bases, size_t nbases);

// Convert from unpacked representation (1 base per byte) to packed
// representation (4 bases per byte)
void binary_seq_pack(uint8_t *restrict ptr,
                     const Nucleotide *restrict bases, size_t len);

// Convert from compact representation (4 bases per byte) to unpacked
// representation (1 base per byte)
void binary_seq_unpack(const uint8_t *restrict ptr,
                       Nucleotide *restrict bases, size_t len);

// Copy a packed path from one place in memory to another, applying left shift
// Shifting by N bases results in N fewer bases in output
// len_bases is length before shifting
// dst needs as many bytes output as input

void binary_seq_cpy_slow(uint8_t *restrict dst, const uint8_t *restrict src,
                         size_t shift, size_t n);

void binary_seq_cpy_med(uint8_t *restrict dst, const uint8_t *restrict src,
                        size_t shift, size_t n);

void binary_seq_cpy_fast(uint8_t *restrict dst, const uint8_t *restrict src,
                         uint8_t shift, size_t n);

#define binary_seq_cpy(dst,src,shift,len) binary_seq_cpy_fast(dst,src,shift,len)

char* binary_seq_to_str(const uint8_t *arr, size_t len, char *str);
void binary_seq_from_str(const char *str, size_t len, uint8_t *arr);

void binary_seq_print(const uint8_t *arr, size_t len, FILE *fout);
void binary_seq_gzprint(const uint8_t *arr, size_t len, gzFile gzout);

// Lexicographic comparison e.g.:
//   a
//   ac
//   ag
//   c
//   cg
//   cgt
// len is in bases
// Returns: -1 if arr0 < arr1, 0 if arr0 == arr1, 1 if arr0 > arr1
int binary_seqs_cmp(const uint8_t *arr0, size_t len0,
                    const uint8_t *arr1, size_t len1);

#endif /* BINARY_SEQ_H_ */
