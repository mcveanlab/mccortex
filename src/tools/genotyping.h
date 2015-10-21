#ifndef GENOTYPING_H_
#define GENOTYPING_H_

#include "binary_kmer.h"
#include "htslib/khash.h"

#include "htslib/vcf.h"

static inline int bk2bits_hash(BinaryKmer bkey) { return binary_kmer_hash(bkey, 0); }
static inline int bk2bits_eq(BinaryKmer k1, BinaryKmer k2) { return binary_kmers_are_equal(k1,k2); }
KHASH_INIT(BkToBits, BinaryKmer, uint64_t, 1, bk2bits_hash, bk2bits_eq);

typedef struct {
  // nkmers are the number of unique kmers that could were in the graph
  // sumcovg are the sum of coverages on those kmers
  // [0] => ref, [1] => alt
  size_t nkmers[2], sumcovg[2];
} VarCovg;

// VCF line
typedef struct
{
  bcf1_t v;
  size_t vidx;
} VcfCovLine;

// Single alt-allele decomposed from a VCF entry
typedef struct {
  VcfCovLine *parent;
  const char *ref, *alt;
  uint32_t pos, reflen, altlen, aid; // variant, allele id
  VarCovg *c; // entry for each colour
  bool has_covg; // if we have fetched coverage into `c`
  size_t nhapk[2]; // number of kmers unique to either allele (0=>ref,1=>alt)
} VcfCovAlt;

typedef struct {
  BinaryKmer bkey;
  uint64_t arbits; // alt-ref-bits
} HaploKmer;

// index after the last base used in this haplotype
#define vcfcovalt_hap_start(v,ks) ((v)->pos <= (ks)-1 ? 0 : (v)->pos - ((ks)-1))
#define vcfcovalt_hap_end(v,ks) ((v)->pos + (v)->reflen + (ks) - 1)

typedef struct GenotyperStruct Genotyper;

#include "madcrowlib/madcrow_buffer.h"
#include "madcrowlib/madcrow_list.h"
madcrow_buffer(haplokmer_buf, HaploKmerBuffer, HaploKmer);
madcrow_list(vc_lines, VcfCovLinePtrList, VcfCovLine*);
madcrow_list(vc_alts,  VcfCovAltPtrList,  VcfCovAlt*);

static inline void vcfcov_alt_wipe_covg(VcfCovAlt *var, size_t ncols)
{
  if(var->has_covg) {
    memset(var->c, 0, ncols*sizeof(var->c[0]));
    memset(var->nhapk, 0, sizeof(var->nhapk));
  }
  var->has_covg = false;
}

int vcfcov_alt_ptr_cmp(const VcfCovAlt *a, const VcfCovAlt *b);
void vcfcov_alts_sort(VcfCovAlt **vars, size_t nvars);

Genotyper* genotyper_init();
void genotyper_destroy(Genotyper *typer);

// Returns non-zero iff kmer occurs in only one of ref/alt in at least one sample
// 0x5 in binary is 0101
#define genotyping_refalt_uniq(b) (((b) ^ ((b)>>1)) & 0x5555555555555555UL)

static inline uint64_t genotyping_refalt_nonuniq(uint64_t b) {
  b = (b & (b>>1)) & 0x5555555555555555UL;
  return b | (b<<1);
}

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
                            const VcfCovAlt *const* vars, size_t nvars,
                            size_t tgtidx, size_t ntgts,
                            const char *chrom, size_t chromlen,
                            size_t kmer_size, HaploKmer **result);

#endif /* GENOTYPING_H_ */
