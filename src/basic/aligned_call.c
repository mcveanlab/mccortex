#include "global.h"
#include "aligned_call.h"

struct CallDecompStruct
{
  nw_aligner_t *nw_aligner;
  scoring_t *scoring;
  alignment_t *aln;
  // TODO: convert output to htslib
  // Output
  FILE *fout;
  DecomposeStats stats;
};

CallDecomp* call_decomp_init(FILE *fout)
{
  CallDecomp *dc = ctx_calloc(1, sizeof(CallDecomp));
  dc->nw_aligner = needleman_wunsch_new();
  dc->aln = alignment_create(1024);
  dc->scoring = ctx_calloc(1, sizeof(dc->scoring[0]));
  scoring_system_default(dc->scoring);
  dc->fout = fout;
  return dc;
}

void call_decomp_destroy(CallDecomp *dc)
{
  alignment_free(dc->aln);
  needleman_wunsch_free(dc->nw_aligner);
  ctx_free(dc->scoring);
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

// Print allele with previous base and no deletions
// 'A--CG-T' with prev_base 'C' => print 'CACGT'
static void print_vcf_allele(int prev_base, const char *allele, size_t len,
                             FILE *fout)
{
  size_t i;
  if(prev_base > 0) fputc((char)prev_base, fout);
  for(i = 0; i < len; i++)
    if(allele[i] != '-')
      fputc(allele[i], fout);
}

// @param vcf_pos is 1-based
// @param prev_base is -1 if SNP otherwise previous base
static void print_vcf_entry(size_t vcf_pos, int prev_base,
                            const char *ref, const char *alt, size_t len,
                            const uint8_t *gts, size_t nsamples,
                            CallDecomp *dc, const AlignedCall *call,
                            size_t max_allele_len)
{
  // Check actual allele length
  size_t i, alt_bases = 0;
  for(i = 0; i < len; i++) alt_bases += (alt[i] != '-');
  if(alt_bases > max_allele_len) { dc->stats.nallele_too_long++; return; }

  // CHROM POS ID REF ALT QUAL FILTER INFO
  fprintf(dc->fout, "%s\t%zu\t.\t", call->chrom->name.b, vcf_pos);
  print_vcf_allele(prev_base, ref, len, dc->fout);
  fputc('\t', dc->fout);
  print_vcf_allele(prev_base, alt, len, dc->fout);
  fputs("\t.\tPASS\t", dc->fout);
  fputs(call->info.b ? call->info.b : ".", dc->fout);
  fputs("\tGT", dc->fout);

  // Print genotypes
  for(i = 0; i < nsamples; i++) {
    fputc('\t', dc->fout);
    fputc(gts[i] ? '1' : '.', dc->fout);
  }

  fputc('\n', dc->fout);

  dc->stats.nvars++;
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

static size_t align_get_len(const char *allele, size_t len)
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
                            const uint8_t *gts, size_t nsamples,
                            CallDecomp *dc, const AlignedCall *call,
                            size_t max_allele_len)
{
  int start, len;
  size_t ref_allele_len, alt_allele_len, ref_pos = call->start;
  int prev_base, vcf_pos;
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

    ref_allele_len = align_get_len(ref, len);
    alt_allele_len = align_get_len(alt, len);
    is_snp = (ref_allele_len == 1 && alt_allele_len == 1);
    vcf_pos = ref_pos+1; // Convert to 1-based

    if(!is_snp) {
      prev_base = ref_pos > 0 ? call->chrom->seq.b[ref_pos-1] : 'N';
      vcf_pos--;
    } else {
      prev_base = -1;
    }

    print_vcf_entry(vcf_pos, prev_base, ref, alt, len, gts, nsamples,
                    dc, call, max_allele_len);

    ref_pos += ref_allele_len;
    ref += len;
    alt += len;
  }
}

void acall_decompose(CallDecomp *dc, const AlignedCall *call,
                     size_t max_allele_len, size_t max_flank_dist)
{
  const char *ref_allele = call->chrom->seq.b + call->start;
  size_t i, ref_len = call->end - call->start;
  const StrBuf *alt;

  ctx_assert2(call->start <= call->end, "%u .. %u", call->start, call->end);

  if(call->end - call->start > max_flank_dist) {
    dc->stats.nflanks_too_far_apart++;
    return; // can't align
  }

  for(i = 0; i < call->n_lines; i++)
  {
    alt = &call->lines[i];
    if(ref_len != alt->end || strncasecmp(ref_allele, alt->b, ref_len))
    {
      needleman_wunsch_align2(ref_allele, alt->b, ref_len, alt->end,
                              dc->scoring, dc->nw_aligner, dc->aln);
      
      align_biallelic(dc->aln->result_a, dc->aln->result_b,
                      call->gts+i*call->n_samples, call->n_samples,
                      dc, call, max_allele_len);
      dc->stats.naligned++;
    }
  }

  // update stats
  dc->stats.ncalls++;
}
