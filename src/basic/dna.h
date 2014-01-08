#ifndef DNA_H_
#define DNA_H_

typedef uint8_t Nucleotide;

extern const char dna_nuc_to_char_arr[4];
extern const Nucleotide dna_char_to_nuc_arr[256];
extern const unsigned char dna_complement_char_arr[256];

#define char_is_acgtn(c) (dna_char_to_nuc_arr[(unsigned char)(c)] <= 4)
#define char_is_acgt(c) (dna_char_to_nuc_arr[(unsigned char)(c)] < 4)

#ifdef NDEBUG
  // Input must be valid n=A,C,G,T,a,c,g,t n=0,1,2,3
  #define dna_nuc_to_char(n)  dna_nuc_to_char_arr[(n)]
  #define dna_char_to_nuc(c)  dna_char_to_nuc_arr[(unsigned char)(c)]
  #define dna_nuc_complement(n) ((Nucleotide)(~(n) & 0x3))
#else
  // Include check for valid input
  #define dna_nuc_to_char(n) ({ assert(((n)&3)==(n)), dna_nuc_to_char_arr[n]; })
  #define dna_char_to_nuc(c)  ({ assert(char_is_acgtn(c)), dna_char_to_nuc_arr[(unsigned char)(c)]; })
  #define dna_nuc_complement(n) ({ assert(((n)&3)==(n)), (Nucleotide)(~(n) & 0x3); })
#endif

#define dna_char_complement(c) (char)dna_complement_char_arr[(unsigned char)(c)]

// length is the length in number of bases
// the char* should have one MORE base than that allocated, to hold '\0'
char* dna_reverse_complement_str(char *str, size_t length);

#endif /* DNA_H_ */
