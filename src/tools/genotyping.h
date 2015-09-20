#ifndef GENOTYPING_H_
#define GENOTYPING_H_

#include "binary_kmer.h"
#include "htslib/khash.h"

static inline int bk2bits_hash(BinaryKmer bkey) { return binary_kmer_hash(bkey, 0); }
static inline int bk2bits_eq(BinaryKmer k1, BinaryKmer k2) { return binary_kmers_are_equal(k1,k2); }
KHASH_INIT(BkToBits, BinaryKmer, uint64_t, 1, bk2bits_hash, bk2bits_eq);

typedef struct {
  const char *ref, *alt;
  size_t pos, reflen, altlen;
} GenoVar;

typedef struct {
  BinaryKmer bkey;
  uint64_t arbits; // alt-ref-bits
} GenoKmer;

#include "madcrowlib/madcrow_buffer.h"
#include "madcrowlib/madcrow_list.h"
madcrow_buffer(genokmer_buf, GenoKmerBuffer, GenoKmer);
madcrow_list(  genovar_list, GenoVarList,  GenoVar);

typedef struct GenotyperStruct Genotyper;

Genotyper* genotyper_init();
void genotyper_destroy(Genotyper *typer);

void genovars_sort(const GenoVar **vars, size_t nvars);

// Returns non-zero iff kmer occurs in only one of ref/alt in at least one sample
// 0x5 in binary is 0101
#define genotyping_refalt_uniq(b) (((b) ^ ((b)>>1)) & 0x5555555555555555UL)

/**
 * Get a list of kmers which support variants.
 * @param typer     initialised memory to use
 * @param vars      variants to genotype and surrounding vars.
 *                  Required sort order: pos, reflen, altlen, alt
 * @param nvars     number of variants in `vars`
 * @param tgtidx    index in vars of first variant to type
 * @param ntgts     number of variants to type
 * @param chrom     reference chromosome
 * @param chromlen  length of reference chromosom
 * @param kmer_size kmer size to type at
 * @return number of kmers
 */
size_t genotyping_get_kmers(Genotyper *typer,
                            const GenoVar **vars, size_t nvars,
                            size_t tgtidx, size_t ntgts,
                            const char *chrom, size_t chromlen,
                            size_t kmer_size, GenoKmer **result);

#endif /* GENOTYPING_H_ */
