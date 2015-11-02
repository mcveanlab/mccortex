#include "global.h"
#include "genotyping.h"
#include "dna.h"
#include "seq_reader.h"

struct GenotyperStruct {
  StrBuf seq;
  khash_t(BkToBits) *h;
  HaploKmerBuffer kmer_buf;
};

#define varend(v) ((v)->pos+(v)->reflen)
#define bits_has(b,i) ((b) & (1UL << (i)))

Genotyper* genotyper_init()
{
  Genotyper *gt = ctx_calloc(1, sizeof(Genotyper));
  strbuf_alloc(&gt->seq, 1024);
  gt->h = kh_init(BkToBits);
  haplokmer_buf_alloc(&gt->kmer_buf, 512);
  return gt;
}

void genotyper_destroy(Genotyper *gt)
{
  strbuf_dealloc(&gt->seq);
  kh_destroy(BkToBits, gt->h);
  haplokmer_buf_dealloc(&gt->kmer_buf);
  ctx_free(gt);
}

int vcfcov_alt_ptr_cmp(const VcfCovAlt *a, const VcfCovAlt *b)
{
  if(a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
  if(a->reflen != b->reflen) return a->reflen < b->reflen ? -1 : 1;
  if(a->altlen != b->altlen) return a->altlen < b->altlen ? -1 : 1;
  return strncmp(a->alt, b->alt, a->altlen);
}

static inline int _vcfcov_alt_ptr_cmp(const void *aa, const void *bb)
{
  const VcfCovAlt *a = *(const VcfCovAlt*const*)aa;
  const VcfCovAlt *b = *(const VcfCovAlt*const*)bb;
  return vcfcov_alt_ptr_cmp(a, b);
}

void vcfcov_alts_sort(VcfCovAlt **vars, size_t nvars)
{
  qsort(vars, nvars, sizeof(*vars), _vcfcov_alt_ptr_cmp);
}

static bool vars_compatible(const VcfCovAlt *const*vars, size_t nvars,
                            uint64_t bits)
{
  ctx_assert(nvars < 64);
  size_t i, end = 0;
  for(i = 0; i < nvars; i++, bits>>=1) {
    if(bits & 1) {
      if(vars[i]->pos < end) return false;
      end = MAX2(end, varend(vars[i]));
    }
  }
  return true;
}

// Generate the DNA string sequence of a haplotype
// Store it in parameter seq
// @param regend not inclusive
static inline void assemble_haplotype_str(StrBuf *seq, const char *chrom,
                                          size_t regstart, size_t regend,
                                          const VcfCovAlt *const*vars,
                                          size_t nvars, uint64_t bits)
{
  strbuf_reset(seq);
  uint64_t i, end = regstart;

  for(i = 0; i < nvars; i++, bits>>=1) {
    if(bits & 1UL) {
      ctx_assert(end <= vars[i]->pos);
      strbuf_append_strn(seq, chrom+end, vars[i]->pos-end);
      strbuf_append_strn(seq, vars[i]->alt, vars[i]->altlen);
      end = vars[i]->pos + vars[i]->reflen;
    }
  }
  strbuf_append_strn(seq, chrom+end, regend-end);
  // printf("hapstr: %s %zu-%zu\n", seq->b, regstart, regend);
}

// Get alt/ref bits from bits specifying which alt alleles are currently used
// We then infer which ref alleles are used and which are missing ref and alt
// We return ref/alt bits which indicate which ref and alt alleles are
// represented
static inline uint64_t altrefbits(const VcfCovAlt *const*vars, size_t nvars,
                                  uint64_t bits)
{
  uint64_t b = 0, hasref;
  size_t i, j, k, vend;

  for(i = j = 0; i < nvars; i++)
  {
    if(bits_has(bits,i)) {
      // set alt bit only
      b |= 1UL << (i*2+1);
      continue;
    }

    // Check if we should set ref bit
    // use j so we don't have to start search from zero each time
    while(j < nvars && varend(vars[j]) <= vars[i]->pos) j++;

    vend = varend(vars[i]);
    hasref = 1;

    for(k = j; k < nvars && vars[k]->pos < vend; k++) {
      if(bits_has(bits,k) && varend(vars[k]) > vars[i]->pos)
      {
        hasref = 0;
        break;
      }
    }
    b |= hasref << (i*2);
  }

  return b;
}

static size_t count_ref_kmers(const char *seq, size_t slen,
                              size_t pos, size_t rlen,
                              size_t ksize)
{
  size_t start = pos < (ksize-1) ? 0 : pos - (ksize-1);
  size_t end = MIN2(pos + rlen + (ksize-1), slen);
  const char *l, *r;
  for(l = seq+pos; l > seq+start && char_is_acgt(*(l-1)); l--) {}
  for(r = seq+pos+rlen; r < seq+end && char_is_acgt(*r); r++) {}
  size_t len = r-l;
  return len < ksize ? 0 : len - ksize + 1;
}

/**
 * Get a list of kmers which support variants.
 *
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
                            const VcfCovAlt *const*vars, size_t nvars,
                            size_t tgtidx, size_t ntgts,
                            const char *chrom, size_t chromlen,
                            size_t kmer_size,
                            HaploKmer **result, uint32_t *nrkmers)
{
  ctx_assert2(0 < nvars && nvars < 64, "nvars: %zu", nvars);
  ctx_assert2(tgtidx < nvars, "tgtidx:%zu >= nvars:%zu ??", tgtidx, nvars);
  ctx_assert2(ntgts <= 32, "Too many targets: %zu", ntgts);

  HaploKmerBuffer *gkbuf = &typer->kmer_buf;
  haplokmer_buf_reset(gkbuf);

  const VcfCovAlt *tgt = vars[tgtidx];

  size_t i, tgtend = tgtidx+ntgts;
  size_t regstart = MIN2(vars[0]->pos, vcfcovalt_hap_start(tgt,kmer_size));
  size_t regend = regstart;

  for(i = 0; i < tgtidx; i++) regend = MAX2(regend, varend(vars[i]));
  for(; i < tgtend; i++) regend = MAX2(regend, vcfcovalt_hap_end(vars[i],kmer_size));
  for(; i < nvars; i++) regend = MAX2(regend, varend(vars[i]));

  regend = MIN2(regend, chromlen);

  StrBuf *seq = &typer->seq;
  khash_t(BkToBits) *h = typer->h;
  BinaryKmer bkmer, bkey;
  int hret;
  khiter_t kiter;

  // Count number of ref kmers
  if(nrkmers) {
    for(i = 0; i < ntgts; i++) {
      nrkmers[i] = count_ref_kmers(chrom+regstart, regend-regstart,
                                   vars[tgtidx+i]->pos-regstart,
                                   vars[tgtidx+i]->reflen,
                                   kmer_size);
    }
  }

  // Faster to clear at the end whilst iterating
  // kh_clear(BkToBits, h);

  // Start with ref haplotype (no variants)
  uint64_t bits, limit, altref_bits;
  size_t cstart, cend, cnext;

  for(bits = 0, limit = 1UL<<nvars; bits < limit; bits++) {
    if(vars_compatible(vars, nvars, bits)) {
      // Construct haplotype
      assemble_haplotype_str(seq, chrom, regstart, regend,
                             vars, nvars, bits);

      altref_bits = altrefbits(vars+tgtidx, ntgts, bits>>tgtidx);

      // Covert to kmers, find/add them to the hash table, OR bits
      // Split contig at non-ACGT bases (e.g. N)
      cnext = 0;
      while((cstart = seq_contig_start2(seq->b, seq->end, NULL, 0,
                                        cnext, kmer_size, 0, 0)) < seq->end)
      {
        cend = seq_contig_end2(seq->b, seq->end, NULL, 0,
                               cstart, kmer_size, 0, 0, &cnext);

        // Get kmers
        bkmer = binary_kmer_from_str(seq->b+cstart, kmer_size);
        bkmer = binary_kmer_right_shift_one_base(bkmer);

        for(i = cstart+kmer_size-1; i < cend; i++) {
          bkmer = binary_kmer_left_shift_add(bkmer, kmer_size,
                                             dna_char_to_nuc(seq->b[i]));
          bkey = binary_kmer_get_key(bkmer, kmer_size);
          kiter = kh_put(BkToBits, h, bkey, &hret);
          if(hret < 0) die("khash table failed: out of memory?");
          if(hret > 0) kh_value(h, kiter) = 0; // initialise if not in table
          kh_value(h, kiter) |= altref_bits;
        }
      }
    }
  }

  size_t nkmers = kh_size(h);
  haplokmer_buf_capacity(gkbuf, nkmers);

  // Fetch and delete every item in the hash table
  for(i = 0, kiter = kh_begin(h); kiter != kh_end(h); ++kiter) {
    if(kh_exist(h, kiter)) {
      bkey = kh_key(h, kiter);
      altref_bits = kh_value(h, kiter);
      kh_del(BkToBits, h, kiter);
      if(genotyping_refalt_uniq(altref_bits)) {
        gkbuf->b[i++] = (HaploKmer){.bkey = bkey, .arbits = altref_bits};
      }
    }
  }
  gkbuf->len = i;

  *result = gkbuf->b;
  return gkbuf->len;
}
