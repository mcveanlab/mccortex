#include "global.h"
#include "genotyping.h"

void genotyper_alloc(Genotyper *typer)
{
  strbuf_alloc(&typer->seq, 1024);
  typer->h = kh_init(BkToBits);
}

void genotyper_dealloc(Genotyper *typer)
{
  strbuf_dealloc(&typer->seq);
  kh_destroy(BkToBits, typer->h);
  memset(typer, 0, sizeof(*typer));
}

int genovar_cmp(const void *aa, const void *bb)
{
  const GenoVar *a = (const GenoVar*)aa, *b = (const GenoVar*)bb;
  if(a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
  if(a->reflen != b->reflen) return a->reflen < b->reflen ? -1 : 1;
  if(a->altlen != b->altlen) return a->altlen < b->altlen ? -1 : 1;
  return strncmp(a->alt, b->alt, a->altlen);
}

void genovars_sort(GenoVar *vars, size_t nvars)
{
  qsort(vars, nvars, sizeof(vars[0]), genovar_cmp);
}

static bool vars_compatible(const GenoVar *vars, size_t nvars, uint64_t bits)
{
  ctx_assert(nvars < 64);
  size_t i, end = 0;
  uint64_t b;
  for(i = 0, b = 1; i < nvars; i++, b<<=1) {
    if(bits & b) {
      if(vars[i].pos < end) return false;
      end = vars[i].pos + vars[i].reflen;
    }
  }
  return true;
}

static inline void assemble_haplotype(StrBuf *seq, const char *chrom,
                                      size_t regstart, size_t regend,
                                      const GenoVar *vars, size_t nvars,
                                      uint64_t bits)
{
  strbuf_reset(seq);
  uint64_t i, end = regstart;

  for(i = 0; i < nvars; i++, bits>>=1) {
    if(bits & 1) {
      ctx_assert(end >= vars[i].pos);
      strbuf_append_strn(seq, chrom+end, vars[i].pos-end);
      strbuf_append_strn(seq, vars[i].alt, vars[i].altlen);
      end = vars[i].pos + vars[i].reflen;
    }
  }
  strbuf_append_strn(seq, chrom+end, regend-end);
}

//                 arararararar r=ref, a=alt
// var:  543210    554433221100
// bits: 010110 -> 011001101001
static inline uint64_t varbits_to_altrefbits(uint64_t bits,
                                             size_t tgtidx,
                                             size_t ntgts)
{
  ctx_assert(ntgts <= 32);
  uint64_t i, r = 0;
  bits >>= tgtidx;
  for(i = 0; i < ntgts; i++, bits>>=1)
    r |= 1UL << (i*2 + (bits&1));
  return r;
}

/**
 * Get coverage on alt allele
 * @param vars must be sorted by pos, then reflen, then altlen, then alt
 * @param colour if -1 population coverage
 */
void genotyping_get_covg(Genotyper *typer,
                         const GenoVar *vars, size_t nvars,
                         size_t tgtidx, size_t ntgts,
                         const char *chrom, size_t chromlen,
                         size_t kmer_size,
                         GenoKmerBuffer *gkbuf)
{
  ctx_assert(nvars < 64);
  ctx_assert(nvars > 0);
  ctx_assert(tgtidx < nvars);
  ctx_assert(ntgts <= 32);

  const GenoVar *tgt = &vars[tgtidx];

  // const size_t kmer_size = db_graph->kmer_size;
  long minpos = MIN2(vars[0].pos, tgt->pos - kmer_size + 1);
  size_t i, regstart, regend;

  regstart = MAX2(minpos, 0);
  regend = regstart;
  for(i = 0; i < nvars; i++) regend = MAX2(regend, vars[i].pos + vars[i].reflen);
  ctx_assert(regend <= chromlen);
  regend = MAX2(regend, tgt->pos + tgt->reflen + kmer_size - 1);
  regend = MIN2(regend, chromlen);

  StrBuf *seq = &typer->seq;
  khash_t(BkToBits) *h = typer->h;
  BinaryKmer bkey;
  int hret;
  khiter_t k;

  kh_clear(BkToBits, h);

  // Start with ref haplotype (no variants)
  uint64_t bits = 0, limit = 1UL<<nvars, altref_bits;

  for(; bits < limit; bits++) {
    if(vars_compatible(vars, nvars, bits)) {
      // Construct haplotype
      assemble_haplotype(seq, chrom, regstart, regend,
                         vars, nvars, bits);

      altref_bits = varbits_to_altrefbits(bits, tgtidx, ntgts);

      // Covert to kmers, find/add them to the hash table, OR bits
      for(i = 0; i + kmer_size <= seq->end; i++) {
        bkey = binary_kmer_from_str(seq->b+i, kmer_size);
        bkey = binary_kmer_get_key(bkey, kmer_size);
        k = kh_put(BkToBits, h, bkey, &hret);
        if(hret < 0) die("khash table failed: out of memory?");
        if(hret > 0) kh_value(h, k) = 0; // initialise if not already in table
        kh_value(h, k) |= altref_bits;
      }
    }
  }

  size_t nkmers = kh_size(h);
  genokmer_buf_capacity(gkbuf, nkmers);

  for(i = 0, k = kh_begin(h); k != kh_end(h); ++k) {
    if(kh_exist(h, k)) {
      bkey = kh_key(h, k);
      altref_bits = kh_value(h, k);
      if(genotyping_refalt_uniq(altref_bits)) {
        gkbuf->b[i++] = (GenoKmer){.bkey = bkey, .arbits = altref_bits};
      }
    }
  }
  gkbuf->len = i;
}
