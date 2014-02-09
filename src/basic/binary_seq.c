#include "global.h"
#include "binary_seq.h"

// byte reverse complement look up table
// Example: byte representing ACTG -> CAGT
//   since A=00, C=01, G=10, T=11
//         ACTG = 00011110 (30)
//   revcmp[30] = 01001011 (0x4B)
//              => CAGT
const uint8_t revcmp_table[256] =
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

void binary_seq_reverse_complement(uint8_t *bases, size_t nbases)
{
  size_t nbytes, top_bits, unused_bits, i, j;
  uint8_t tmpbases;

  nbytes = (nbases+3)/4;
  top_bits = 2*bases_in_top_byte(nbases);
  unused_bits = 8 - top_bits;

  if(nbases == 0) return;
  if(nbytes == 1) { bases[0] = revcmp_table[bases[0]] >> unused_bits; return; }

  for(i = 0, j = nbytes-1; i <= j; i++, j--) {
    tmpbases = bases[i];
    bases[i] = revcmp_table[bases[j]];
    bases[j] = revcmp_table[tmpbases];
  }

  // shift
  if(unused_bits > 0) {
    for(i = 0; i+1 < nbytes; i++)
      bases[i] = (bases[i]>>unused_bits) | (bases[i+1]<<top_bits);
    bases[nbytes-1] >>= unused_bits;
  }
}

// Convert from unpacked representation (1 base per byte) to packed
// representation (4 bases per byte)
void binary_seq_pack(uint8_t *restrict ptr,
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
void binary_seq_unpack(const uint8_t *restrict ptr,
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

void binary_seq_cpy_slow(uint8_t *restrict dst, const uint8_t *restrict src,
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

void binary_seq_cpy_med(uint8_t *restrict dst, const uint8_t *restrict src,
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
void binary_seq_cpy_fast(uint8_t *restrict dst, const uint8_t *restrict src,
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

void binary_seq_to_str(const uint8_t *arr, size_t len, char *str)
{
  const char *end = str+len;
  size_t b = 0, o = 0;
  for(; str < end; str++, o+=2, b += (o == 8), o &= 7) {
    *str = dna_nuc_to_char((arr[b] >> o) & 3);
  }
  *str = '\0';
}
