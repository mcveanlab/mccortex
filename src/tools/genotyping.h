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

  // Covg
  size_t refkmers, refsumcovg;
  size_t altkmers, altsumcovg;
} GenoVar;

typedef struct {
  BinaryKmer bkey;
  uint64_t arbits; // alt-ref-bits
} GenoKmer;

#include "madcrowlib/madcrow_buffer.h"
#include "madcrowlib/madcrow_list.h"
madcrow_buffer(genokmer_buf, GenoKmerBuffer, GenoKmer);
// madcrow_buffer(genovar_buf,  GenoVarBuffer,  GenoVar);
madcrow_list(  genovar_list, GenoVarList,  GenoVar);

typedef struct {
  StrBuf seq;
  khash_t(BkToBits) *h;
  GenoKmerBuffer kmer_buf;
} Genotyper;

void genotyper_alloc(Genotyper *typer);
void genotyper_dealloc(Genotyper *typer);

void genovars_sort(GenoVar *vars, size_t nvars);

// Returns 1 if kmer occurs in ref/alt only for a given haplotype
// print out if pairs of bits are 01 or 10, not 11, 00
// 0x5 in binary is 0101
#define genotyping_refalt_uniq(b) (((b) ^ ((b)>>1)) & 0x5555555555555555UL)

/**
 * Get coverage on alt allele. Results stored in typer->kmer_buf.
 * @param vars must be sorted by pos, then reflen, then altlen, then alt
 * @param colour if -1 population coverage
 */
void genotyping_get_covg(Genotyper *typer,
                         const GenoVar *vars, size_t nvars,
                         size_t tgtidx, size_t ntgts,
                         const char *chrom, size_t chromlen,
                         size_t kmer_size);

#endif /* GENOTYPING_H_ */
