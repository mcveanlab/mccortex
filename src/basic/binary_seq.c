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
  0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f,
  0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0x0f,
  0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b,
  0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0x0b,
  0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27,
  0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x07,
  0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23,
  0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x03,
  0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e,
  0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0x0e,
  0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a,
  0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0x0a,
  0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26,
  0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x06,
  0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22,
  0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x02,
  0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d,
  0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0x0d,
  0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29,
  0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x09,
  0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25,
  0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x05,
  0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21,
  0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x01,
  0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c,
  0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0x0c,
  0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28,
  0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x08,
  0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24,
  0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x04,
  0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20,
  0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x00
};

// Leave top bases intact: e.g. bases=TACG revcmp(bases,2): ATCG
void binary_seq_reverse_complement(uint8_t *bases, size_t nbases)
{
  if(!nbases) return;

  size_t nbytes, top_bits, unused_bits, i, j;
  uint8_t tmpbases, topbyte;

  nbytes = binary_seq_mem(nbases);
  top_bits = 2*bases_in_top_byte(nbases);
  unused_bits = 8 - top_bits;
  topbyte = bases[nbytes-1];

  if(nbytes == 1) {
    // Cast to size_t so >>8 isn't undefined
    bases[0] = ((size_t)revcmp_table[bases[0]] >> unused_bits) |
               (topbyte & (((size_t)255 >> top_bits) << top_bits));
    return;
  }

  for(i = 0, j = nbytes-1; i <= j; i++, j--) {
    tmpbases = bases[i];
    bases[i] = revcmp_table[bases[j]];
    bases[j] = revcmp_table[tmpbases];
  }

  // shift
  if(unused_bits > 0) {
    for(i = 0; i+1 < nbytes; i++)
      bases[i] = (bases[i]>>unused_bits) | (bases[i+1]<<top_bits);
    bases[nbytes-1] = (bases[nbytes-1] >> unused_bits) |
                      (topbyte & ((255>>top_bits)<<top_bits));
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
  m = binary_seq_mem(n);
  dstn = binary_seq_mem(n-shift);
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
  m = binary_seq_mem(n);
  dstn = binary_seq_mem(n-shift);
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
  size_t src_bytes = binary_seq_mem(n), dst_bytes = binary_seq_mem(n-shift);

  if(shift >= n) { dst[0] = 0; return; }
  if(!shift) {
    memcpy(dst, src, src_bytes);
    dst[src_bytes-1] &= bitmask64(n*2-(src_bytes-1)*8); // mask top byte
    return;
  }

  size_t nwords64, endbyte, byte, bitshift;
  uint64_t word = 0;

  bitshift = shift*2;
  nwords64 = (n-1)/32; // -1 so we can look ahead
  endbyte = nwords64*8;
  dst[dst_bytes-1] = 0; // zero top byte, so we don't complain when masking later

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

char* binary_seq_to_str(const uint8_t *arr, size_t len, char *str)
{
  char *ptr = str;
  const char *end = str+len;
  size_t b = 0, o = 0;
  while(ptr < end) {
    *ptr = dna_nuc_to_char((arr[b] >> o) & 3);
    ptr++; o += 2; b += (o == 8); o &= 7;
  }
  *ptr = '\0';
  return str;
}

void binary_seq_from_str(const char *str, size_t len, uint8_t *arr)
{
  const char *end = str+len;
  size_t b = 0, o = 0;
  memset(arr, 0, binary_seq_mem(len));
  while(str < end) {
    arr[b] |= dna_char_to_nuc(*str) << o;
    str++; o += 2; b += (o == 8); o &= 7;
  }
}

void binary_seq_print(const uint8_t *arr, size_t len, FILE *fout)
{
  size_t i, b = 0, o = 0;
  for(i = 0; i < len; i++) {
    fputc(dna_nuc_to_char((arr[b] >> o) & 3), fout);
    o+=2; b += (o == 8); o &= 7;
  }
}

void binary_seq_gzprint(const uint8_t *arr, size_t len, gzFile gzout)
{
  size_t i, b = 0, o = 0;
  for(i = 0; i < len; i++) {
    gzputc(gzout, dna_nuc_to_char((arr[b] >> o) & 3));
    o+=2; b += (o == 8); o &= 7;
  }
}

// Lexicographic comparison
// len is in bases
// Returns: -1 if arr0 < arr1, 0 if arr0 == arr1, 1 if arr0 > arr1
int binary_seqs_cmp(const uint8_t *arr0, size_t len0,
                    const uint8_t *arr1, size_t len1)
{
  size_t i, b = 0, o = 0, len = MIN2(len0, len1);
  int x, y, ret;
  // Loop over one base at a time
  for(i = 0; i < len; i++) {
    x = (arr0[b] >> o) & 3;
    y = (arr1[b] >> o) & 3;
    if((ret = x - y) != 0) return ret;
    o+=2; b += (o == 8); o &= 7;
  }
  return (long)len0 - (long)len1;
}
