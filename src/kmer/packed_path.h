#ifndef PACKED_PATH_H_
#define PACKED_PATH_H_

#include "cortex_types.h"
#include "dna.h"
#include "bit_macros.h"

typedef uint64_t PathIndex;
typedef uint16_t PathLen;

#define PATH_NULL UINT64_MAX
#define PATH_LEN_BITS 15
#define MAX_PATHLEN ((1UL<<PATH_LEN_BITS)-1)

//
// Functions for PackedPaths
//
// PackedPath doesn't have a struct, instead is laid out as:
//   {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}
// Where N=round_up(num_of_colours/8), M=round_up(len/4)
// len holds both the length (lower 15 bits), and orientations (top bit)

// Number of bytes a packed path takes up
#define packedpath_mem2(colbytes,nbytes) \
        (sizeof(PathIndex)+(colbytes)+sizeof(PathLen)+(nbytes))

#define packedpath_mem(p,colbytes) \
        packedpath_mem2(colbytes,packedpath_pbytes(p,colbytes))

#define packedpath_pbytes(p,colbytes) \
        packedpath_len_nbytes(packedpath_get_len(p,colbytes))

//
// Previous
//

static inline PathIndex packedpath_get_prev(const uint8_t *ptr) {
  PathIndex idx; memcpy(&idx, ptr, sizeof(PathIndex)); return idx;
}

static inline void packedpath_set_prev(uint8_t *ptr, PathIndex idx) {
  memcpy(ptr, &idx, sizeof(PathIndex));
}

//
// Colour bitset
//
#define packedpath_del_col(p,col) bitset_del((p)+sizeof(PathIndex),col)
#define packedpath_set_col(p,col) bitset_set((p)+sizeof(PathIndex),col)
#define packedpath_has_col(p,col) bitset_get((p)+sizeof(PathIndex),col)

// Get pointer to colset
#define packedpath_get_colset(ptr) ((ptr) + sizeof(PathIndex))

//
// Length and Orientation
//

// Number of bytes needed to store n bases (2 bits per base, 4 per byte)
#define packedpath_len_nbytes(nbases) (((nbases)+3)/4)

#define PP_ORIENTMASK (1U << PATH_LEN_BITS)
#define PP_LENMASK (PP_ORIENTMASK-1)

#define packedpath_len(w) ((w) & PP_LENMASK)
#define packedpath_or(w) ((w) >> PATH_LEN_BITS)

#define packedpath_combine_lenorient(len,orient) \
        ((PathLen)(((PathLen)(orient))<<PATH_LEN_BITS)|(len))

static inline PathLen packedpath_get_lenword(const uint8_t *ptr, size_t colbytes)
{
  PathLen len; memcpy(&len, ptr+sizeof(PathIndex)+colbytes, sizeof(PathLen));
  return len;
}

static inline void packedpath_set_lenword(uint8_t *ptr, size_t colbytes,
                                          PathLen len_orient)
{
  memcpy(ptr+sizeof(PathIndex)+colbytes, &len_orient, sizeof(PathLen));
}

static inline PathLen packedpath_get_len(const uint8_t *ptr, size_t colbytes) {
  return packedpath_len(packedpath_get_lenword(ptr,colbytes));
}

static inline Orientation packedpath_get_orient(const uint8_t *ptr, size_t colbytes) {
  return packedpath_or(packedpath_get_lenword(ptr,colbytes));
}

static inline void packedpath_set_len(uint8_t *ptr, size_t colbytes, PathLen len)
{
  // Combine with orientation
  len |= packedpath_get_lenword(ptr, colbytes) & PP_ORIENTMASK;
  packedpath_set_lenword(ptr, colbytes, len);
}

static inline void packedpath_set_orient(uint8_t *ptr, size_t colbytes,
                                         Orientation orient)
{
  // Combine with orientation
  PathLen len, len_orient;
  len = packedpath_get_lenword(ptr, colbytes) & PP_LENMASK;
  len_orient = packedpath_combine_lenorient(len, orient);
  packedpath_set_lenword(ptr, colbytes, len_orient);
}

// Get length and orientation together
static inline PathLen packedpath_get_len_orient(const uint8_t *ptr,
                                                size_t colbytes, PathLen *len,
                                                Orientation *orient)
{
  PathLen lenword = packedpath_get_lenword(ptr, colbytes);
  *len = packedpath_len(lenword);
  *orient = packedpath_or(lenword);
  return lenword;
}

static inline void packedpath_set_len_orient(uint8_t *ptr, size_t colbytes,
                                             PathLen len, Orientation orient)
{
  PathLen len_orient = packedpath_combine_lenorient(len, orient);
  packedpath_set_lenword(ptr, colbytes, len_orient);
}

//
// Path
//

// Get pointer to first byte of path
#define packedpath_seq(p,colbytes) ((p)+sizeof(PathIndex)+colbytes+sizeof(PathLen))


// Convert from unpacked representation (1 base per byte) to packed
// representation (4 bases per byte)
static inline void pack_bases(uint8_t *restrict ptr,
                              const Nucleotide *restrict bases, size_t len)
{
  size_t i, full_bytes = len/4;
  const uint8_t *endptr = ptr+full_bytes;

  for(i = 0; ptr < endptr; ptr++, i += 4) {
    *ptr = bases[i] | (uint8_t)(bases[i+1]<<2) |
           (uint8_t)(bases[i+2]<<4) | (uint8_t)(bases[i+3]<<6);
  }

  // Do last byte
  if(len & 3) {
    *ptr = 0;
    switch(len & 3) {
      case 3: *ptr = bases[--len];
      case 2: *ptr = (uint8_t)((*ptr)<<2) | bases[--len];
      case 1: *ptr = (uint8_t)((*ptr)<<2) | bases[--len];
    }
  }
}

// Convert from compact representation (4 bases per byte) to unpacked
// representation (1 base per byte)
static inline void unpack_bases(const uint8_t *restrict ptr,
                                Nucleotide *restrict bases, size_t len)
{
  size_t i, full_bytes = len/4;
  const uint8_t *endptr = ptr+full_bytes;

  for(i = 0; ptr < endptr; ptr++) {
    bases[i++] =  (*ptr)     & 3;
    bases[i++] = ((*ptr)>>2) & 3;
    bases[i++] = ((*ptr)>>4) & 3;
    bases[i++] = ((*ptr)>>6);
  }

  // Do last byte
  switch(len & 3) {
    case 3: bases[--len] = ((*ptr)>>4) & 3;
    case 2: bases[--len] = ((*ptr)>>2) & 3;
    case 1: bases[--len] = ((*ptr))    & 3;
  }
}

// Copy a packed path from one place in memory to another, applying left shift
// Shifting by N bases results in N fewer bases in output
// len_bases is length before shifting
// dst needs as many bytes output as input

static inline void packed_cpy_slow(uint8_t *restrict dst,
                                   const uint8_t *restrict src,
                                   size_t shift, size_t n)
{
  size_t i, m, sb, dstn;
  if(shift >= n) { dst[0] = 0; return; }
  sb = shift*2;
  m = (n+3)/4;
  dstn = (n-shift+3)/4;
  dst[dstn-1] = 0;
  for(i = 0; i+1 < m; i++) dst[i] = (src[i]>>sb) | (uint8_t)(src[i+1]<<(8-sb));
  dst[dstn-1] |= src[dstn-1] >> sb;
  dst[dstn-1] &= bitmask64((n-shift)*2-(dstn-1)*8); // mask top byte
}

static inline void packed_cpy_med(uint8_t *restrict dst,
                                  const uint8_t *restrict src,
                                  size_t shift, size_t n)
{
  size_t i, m, sb, dstn;
  if(shift >= n) { dst[0] = 0; return; }
  sb = shift*2;
  m = (n+3)/4;
  dstn = (n-shift+3)/4;
  dst[dstn-1] = 0;
  switch(shift) {
    case 0: memcpy(dst, src, m); break;
    case 1: for(i=0;i+1<m;i++){ dst[i] = (src[i]>>2) | (uint8_t)(src[i+1]<<6); } break;
    case 2: for(i=0;i+1<m;i++){ dst[i] = (src[i]>>4) | (uint8_t)(src[i+1]<<4); } break;
    case 3: for(i=0;i+1<m;i++){ dst[i] = (src[i]>>6) | (uint8_t)(src[i+1]<<2); } break;
  }
  dst[dstn-1] |= src[dstn-1] >> sb;
  dst[dstn-1] &= bitmask64((n-shift)*2-(dstn-1)*8); // mask top byte
}

// Copy 8 bytes at a time
static inline void packed_cpy_fast(uint8_t *restrict dst,
                                   const uint8_t *restrict src,
                                   uint8_t shift, size_t n)
{
  size_t src_bytes = (n+3)/4, dst_bytes = (n-shift+3)/4;

  if(shift >= n) { dst[0] = 0; return; }
  if(!shift) {
    memcpy(dst, src, src_bytes);
    dst[src_bytes-1] &= bitmask64(n*2-(src_bytes-1)*8); // mask top byte
    return;
  }

  size_t nwords64, endbyte, byte, bitshift;
  uint64_t word;

  bitshift = shift*2;
  nwords64 = (n-1)/32; // -1 so we can look ahead
  endbyte = nwords64*8;

  for(byte=0; byte<endbyte; byte+=8) {
    memcpy(&word, &src[byte], 8);
    word = (word >> bitshift) | ((uint64_t)src[byte+8] << (64-bitshift));
    memcpy(&dst[byte], &word, 8);
  }

  if(byte < dst_bytes) {
    size_t rem_src_bytes = src_bytes - byte;
    size_t rem_dst_bytes = dst_bytes - byte;

    memcpy(&word, &src[byte], rem_src_bytes);
    word >>= bitshift;
    memcpy(&dst[byte], &word, rem_dst_bytes);
  }

  dst[dst_bytes-1] &= bitmask64((n-shift)*2-(dst_bytes-1)*8);
}

#define packed_cpy(dst,src,shift,len) packed_cpy_fast(dst,src,shift,len)

// Fetch a given base. Four bases per byte
// ptr should point directly to sequence
static inline Nucleotide packed_fetch(const uint8_t *ptr, size_t idx)
{
  size_t byte = idx / 4, offset = (idx & 3)*2;
  return (ptr[byte] >> offset) & 3;
}

// Add a given base. Four bases per byte
// ptr should point directly to sequence
// Doesn't provide any masking - should already be zeroed
static inline void packed_add(uint8_t *ptr, size_t idx, Nucleotide nuc)
{
  size_t byte = idx / 4, offset = (idx & 3)*2;
  ptr[byte] |= nuc << offset;
}

#endif /* PACKED_PATH_H_ */
