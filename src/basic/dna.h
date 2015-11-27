#ifndef DNA_H_
#define DNA_H_

typedef uint8_t Nucleotide;

#define STRAND_PLUS 0
#define STRAND_MINUS 1

extern const char dna_nuc_to_char_arr[4];
extern const Nucleotide dna_char_to_nuc_arr[256];
extern const unsigned char dna_complement_char_arr[256];

extern const char dna_char_to_vcf_char[256]; // all non-ACGT chars are 'N'

#define char_is_acgtn(c) (dna_char_to_nuc_arr[(uint8_t)(c)] <= 4)
#define char_is_acgt(c) (dna_char_to_nuc_arr[(uint8_t)(c)] < 4)

#define char_to_vcf_char(c) (dna_char_to_vcf_char[(uint8_t)(c)])

// Include check for valid input
#define dna_nuc_to_char(n) ({ ctx_assert(((n)&3)==(n)); dna_nuc_to_char_arr[n]; })
#define dna_nuc_complement(n) ({ ctx_assert(((n)&3)==(n)); (Nucleotide)(~(n) & 0x3); })
#define dna_char_to_nuc(c)  ({ ctx_assert2(char_is_acgt(c),"%c",c); dna_char_to_nuc_arr[(uint8_t)(c)]; })
#define dna_char_complement(c) ({ ctx_assert2(char_is_acgt(c),"%i",c); (char)dna_complement_char_arr[(uint8_t)(c)]; })

#define dna_reverse_complement_str(str,len) dna_revcomp_str(str,str,len)

/**
 * Reverse complement a string, copying the result into a different memory
 * location. src,dst can point to the same string.
 * DOES NOT NULL TERMINATE dst
 *
 * @param length is the length in number of bases to reverse complement
 * @param dst read characters from dst
 * @param src result written to memory pointed to by src
 * @return pointer to dst
**/
char* dna_revcomp_str(char *dst, const char *src, size_t length);

// Generate a random dna str "ACGT" of length `len`, terminated with a \0 at
// position `len`. `str` must be at least of size `len`+1.
// Useful for testing
char* dna_rand_str(char *str, size_t len);

// compare a with the reverse complement of b
int dna_revncasecmp(const char *a, const char *b, size_t len);

// out must be at least 11 bytes long: "A, C, G, T"
size_t dna_bases_list_to_str(const bool bases[4], char *out);

// Case insensitive comparison that converts non-ACGT characters to N before
// comparing. Useful for VCF ref comparisons.
// Returns true iff sequences match
static inline bool dna_ref_vcf_match(const char *chr, const char *ref, size_t len)
{
  const char *end;
  for(end = chr+len; chr < end; chr++, ref++) {
    ctx_assert2(*ref, "Ref too short");
    if(char_to_vcf_char(*chr) != char_to_vcf_char(*ref))
      return false;
  }
  return true;
}

#endif /* DNA_H_ */
