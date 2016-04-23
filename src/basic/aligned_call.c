#include "global.h"
#include "aligned_call.h"
#include "dna.h"
#include "carrays/carrays.h"

//
// AlignedCall
//

void acall_resize(AlignedCall *call, size_t n_lines, size_t n_samples)
{
  gca_resize(call->lines, call->s_lines, n_lines);
  gca_resize(call->gts, call->s_gts, n_lines*n_samples);
  call->n_lines = n_lines;
  call->n_samples = n_samples;
  call->n_gts = n_lines*n_samples;
}

void acall_destroy(AlignedCall *call)
{
  size_t i;
  for(i = 0; i < call->n_lines; i++) strbuf_dealloc(&call->lines[i]);
  free(call->lines);
  free(call->gts);
  strbuf_dealloc(&call->info);
  ctx_free(call);
}


//
// Decomposer
//

struct CallDecompStruct
{
  nw_aligner_t *nw_aligner;
  scoring_t *scoring;
  alignment_t *aln;
  htsFile *vcffh;
  bcf_hdr_t *vcfhdr;
  bcf1_t *v;
  StrBuf sbuf;
  DecomposeStats stats;
};

CallDecomp* call_decomp_init(htsFile *vcffh, bcf_hdr_t *vcfhdr)
{
  CallDecomp *dc = ctx_calloc(1, sizeof(CallDecomp));
  dc->nw_aligner = needleman_wunsch_new();
  dc->aln = alignment_create(1024);
  dc->scoring = ctx_calloc(1, sizeof(dc->scoring[0]));
  scoring_system_default(dc->scoring);
  dc->vcffh = vcffh;
  dc->vcfhdr = vcfhdr;
  dc->v = bcf_init();
  strbuf_alloc(&dc->sbuf, 256);
  return dc;
}

void call_decomp_destroy(CallDecomp *dc)
{
  alignment_free(dc->aln);
  needleman_wunsch_free(dc->nw_aligner);
  ctx_free(dc->scoring);
  bcf_destroy(dc->v);
  strbuf_dealloc(&dc->sbuf);
  ctx_free(dc);
}

scoring_t* call_decomp_get_scoring(CallDecomp *dc)
{
  return dc->scoring;
}

void call_decomp_cpy_stats(DecomposeStats *stats, const CallDecomp *dc)
{
  memcpy(stats, &dc->stats, sizeof(*stats));
}

//
// Decompose AlignedCall
//

// Strip indels ('-') from allele and add to string buffer
static inline void print_vcf_allele(const char *allele, size_t len,
                                    int8_t prev_base, int8_t next_base,
                                    StrBuf *sbuf)
{
  size_t i;
  if(prev_base > 0) strbuf_append_char(sbuf, char_to_vcf_char(prev_base));
  strbuf_ensure_capacity(sbuf, sbuf->end+len);
  for(i = 0; i < len; i++) {
    if(allele[i] != '-')
      sbuf->b[sbuf->end++] = char_to_vcf_char(allele[i]);
  }
  sbuf->b[sbuf->end] = 0;
  if(next_base > 0) strbuf_append_char(sbuf, char_to_vcf_char(next_base));
}

// @param vcf_pos is 0-based
// @param prev_base is -1 if SNP otherwise previous base
// @param next_base is -1 unless indel at position 0
static void print_vcf_entry(size_t vcf_pos, int8_t prev_base, int8_t next_base,
                            const char *ref, const char *alt, size_t len,
                            const uint8_t *gts, size_t nsamples,
                            CallDecomp *dc, const AlignedCall *call,
                            size_t max_allele_len)
{
  dc->stats.nvars++;

  StrBuf *sbuf = &dc->sbuf;
  strbuf_reset(sbuf);

  // Check actual allele length
  size_t i, alt_bases = 0;
  for(i = 0; i < len; i++) alt_bases += (alt[i] != '-');
  if(alt_bases > max_allele_len) { dc->stats.nallele_too_long++; return; }

  // CHROM POS ID REF ALT QUAL FILTER INFO
  strbuf_append_str(sbuf, call->chrom->name.b);
  strbuf_append_char(sbuf, '\t');
  strbuf_append_ulong(sbuf, vcf_pos+1);
  strbuf_append_str(sbuf, "\t.\t");
  print_vcf_allele(ref, len, prev_base, next_base, sbuf);
  strbuf_append_char(sbuf, '\t');
  print_vcf_allele(alt, len, prev_base, next_base, sbuf);
  strbuf_append_str(sbuf, "\t.\tPASS\t");
  strbuf_append_str(sbuf, call->info.b ? call->info.b : ".");
  strbuf_append_str(sbuf, "\tGT");

  // Print genotypes
  for(i = 0; i < nsamples; i++) {
    strbuf_append_char(sbuf, '\t');
    strbuf_append_char(sbuf, gts[i] ? '1' : '.');
  }

  // fprintf(stderr, " prev_base:%i next_base:%i info:%s\n", prev_base, next_base, call->info.b);
  // fprintf(stderr, "%s [%zu vs %zu]\n", sbuf->b, sbuf->end, strlen(sbuf->b));

  kstring_t ks = {.l = sbuf->end, .m = sbuf->size, .s = sbuf->b};
  if(vcf_parse(&ks, dc->vcfhdr, dc->v) != 0)
    die("Cannot construct VCF entry: %s", sbuf->b);
  if(bcf_write(dc->vcffh, dc->vcfhdr, dc->v) != 0)
    die("Cannot write VCF entry [nsamples: %zu vs %zu]", nsamples, (size_t)bcf_hdr_nsamples(dc->vcfhdr));
  // Move back into our string buffer
  sbuf->b = ks.s;
  sbuf->size = ks.m;

  dc->stats.nvars_printed++;
}

// `ref` and `alt` are aligned alleles - should both be same length strings
// of 'ACGT-'
// return first mismatch position or -1
static int align_get_start(const char *ref, const char *alt)
{
  const char *start = ref;
  while(*ref) {
    if(*ref != *alt) return (ref - start);
    ref++; alt++;
  }
  return -1;
}

// `ref` and `alt` are aligned alleles - should both be same length strings
// of 'ACGT-'
// return first matching position
static int align_get_end(const char *ref, const char *alt)
{
  int i = 0;
  while(ref[i] && ref[i] != alt[i]) i++;
  return i;
}

static size_t align_get_nbases(const char *allele, size_t len)
{
  size_t i, nbases = 0;
  for(i = 0; i < len; i++)
    if(allele[i] != '-')
      nbases++;
  return nbases;
}

/**
 * @param ref_pos is 0-based here
 * @param info is extra text to print in the info field of each variant (may be NULL)
 * @param genotypes is strings to print in genotypes columns, of length num_samples.
 *                  It may be NULL.
 * @return number of variants printed
 */
static void align_biallelic(const char *ref, const char *alt,
                            const read_t *chrom,
                            const uint8_t *gts, size_t nsamples,
                            CallDecomp *dc, const AlignedCall *call,
                            size_t max_allele_len)
{
  int32_t start, len;
  size_t ref_nbases, alt_nbases, ref_pos = call->start, ref_end, vcf_pos;
  int8_t prev_base, next_base;
  bool is_snp;

  // printf("--\n ref: %s\n alt: %s\n", ref, alt);

  while((start = align_get_start(ref, alt)) > -1)
  {
    ref_pos += start; // assume ref[i]==alt[i] means ref[i]!='-'
    ref += start;
    alt += start;
    len = align_get_end(ref, alt);

    // printf("ref: %.*s\nalt: %.*s\nref_pos: %zu start: %i len %i\n",
    //        len, ref, len, alt, ref_pos, start, len);

    ref_nbases = align_get_nbases(ref, len);
    alt_nbases = align_get_nbases(alt, len);
    is_snp = (ref_nbases == 1 && alt_nbases == 1);
    ref_end = ref_pos+ref_nbases;
    vcf_pos = ref_pos; // copy in case we need left padding base

    // If one allele is going to be empty, we need a padding base
    // If ref_pos == 0, add extra base to end instead
    prev_base = next_base = -1;
    if(!is_snp) {
      if(ref_pos > 0) prev_base = chrom->seq.b[--vcf_pos];
      else if(ref_end < chrom->seq.end) next_base = chrom->seq.b[ref_end];
    }

    if(is_snp || prev_base > 0 || next_base > 0) {
      print_vcf_entry(vcf_pos, prev_base, next_base, ref, alt, len,
                      gts, nsamples, dc, call, max_allele_len);
    }

    ref_pos += ref_nbases;
    ref += len;
    alt += len;
  }
}

void acall_decompose(CallDecomp *dc, const AlignedCall *call,
                     size_t max_line_len, size_t max_allele_len)
{
  dc->stats.ncalls++;
  if(call->chrom == NULL) { return; }
  dc->stats.ncalls_mapped++;

  const read_t *chrom = call->chrom;
  const char *ref_allele = chrom->seq.b + call->start;
  size_t i, ref_len = call->end - call->start;
  const StrBuf *alt;

  ctx_assert2(call->start <= call->end, "%u .. %u", call->start, call->end);

  if(ref_len > max_line_len) {
    dc->stats.ncalls_ref_allele_too_long++;
    return; // can't align
  }

  dc->stats.nlines += call->n_lines;

  // printf("chr:%s %u - %u\n", call->chrom->name.b, call->start, call->end);

  for(i = 0; i < call->n_lines; i++)
  {
    alt = &call->lines[i];
    ctx_assert(strlen(alt->b) == alt->end);

    // Quick check if sequence too long or are matching
    if(alt->end > max_line_len) {
      dc->stats.nlines_too_long++;
    } else if(ref_len == alt->end && strncasecmp(ref_allele, alt->b, ref_len) == 0) {
      dc->stats.nlines_match_ref++;
    } else {
      // printf("REF: '%*.s' [%zu]\n", (int)ref_len, ref_allele, ref_len);
      // printf("ALT: '%*.s' [%zu]\n", (int)alt->end, alt->b, alt->end);

      needleman_wunsch_align2(ref_allele, alt->b, ref_len, alt->end,
                              dc->scoring, dc->nw_aligner, dc->aln);

      // printf("ALNA: %s\n", dc->aln->result_a);
      // printf("ALNB: %s\n", dc->aln->result_b);

      align_biallelic(dc->aln->result_a, dc->aln->result_b, chrom,
                      call->gts+i*call->n_samples, call->n_samples,
                      dc, call, max_allele_len);
      dc->stats.nlines_mapped++;
    }
  }
}
