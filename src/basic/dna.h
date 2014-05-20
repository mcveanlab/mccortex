#ifndef DNA_H_
#define DNA_H_

typedef uint8_t Nucleotide;

#define STRAND_PLUS 0
#define STRAND_MINUS 1

extern const char dna_nuc_to_char_arr[4];
extern const Nucleotide dna_char_to_nuc_arr[256];
extern const unsigned char dna_complement_char_arr[256];

#define char_is_acgtn(c) (dna_char_to_nuc_arr[(unsigned char)(c)] <= 4)
#define char_is_acgt(c) (dna_char_to_nuc_arr[(unsigned char)(c)] < 4)

// Include check for valid input
#define dna_nuc_to_char(n) ({ ctx_assert(((n)&3)==(n)), dna_nuc_to_char_arr[n]; })
#define dna_nuc_complement(n) ({ ctx_assert(((n)&3)==(n)), (Nucleotide)(~(n) & 0x3); })
#define dna_char_to_nuc(c)  ({ ctx_assert2(char_is_acgt(c),"%c",c), dna_char_to_nuc_arr[(unsigned char)(c)]; })
#define dna_char_complement(c) ({ ctx_assert2(char_is_acgt(c),"%c",c), (char)dna_complement_char_arr[(unsigned char)(c)]; })

// length is the length in number of bases
// the char* should have one MORE base than that allocated, to hold '\0'
char* dna_reverse_complement_str(char *str, size_t length);

// Generate a random dna str "ACGT" of length cap-1, terminated with a \0 at
// position cap-1. If cap is 0, does nothing. Useful for testing
char* dna_rand_str(char *str, size_t cap);

#endif /* DNA_H_ */
