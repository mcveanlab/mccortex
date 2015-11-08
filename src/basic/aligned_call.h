#ifndef ALIGNED_CALL_
#define ALIGNED_CALL_

#include "carrays/carrays.h"
#include "seq-align/src/needleman_wunsch.h"
#include "seq_file/seq_file.h"


// Variant haplotypes aligned to a region of the reference
// if not aligned, chrom is NULL, other fields are undefined
typedef struct
{
  StrBuf *lines;
  size_t n_lines, s_lines; // s_lines is size of lines array
  const read_t *chrom; // set to NULL if not aligned
  uint32_t start, end;

  StrBuf info;
  size_t n_samples, s_gts, n_gts;
  uint8_t *gts; // gts[line*n_gts + sample] = 0 (unknown) or 1 (ref)
} AlignedCall;


#define acall_init() ctx_calloc(1, sizeof(AlignedCall))

static inline void acall_resize(AlignedCall *call, size_t n_lines, size_t n_samples)
{
  gca_resize(call->lines, call->s_lines, n_lines);
  gca_resize(call->gts, call->s_gts, n_lines*n_samples);
  call->n_lines = n_lines;
  call->n_samples = n_samples;
  call->n_gts = n_lines*n_samples;
}

static inline void acall_destroy(AlignedCall *call)
{
  size_t i;
  for(i = 0; i < call->n_lines; i++) strbuf_dealloc(&call->lines[i]);
  free(call->lines);
  free(call->gts);
  strbuf_dealloc(&call->info);
  ctx_free(call);
}

//
// Decompose aligned calls to make a VCF
//
typedef struct CallDecompStruct CallDecomp;

typedef struct {
  uint64_t ncalls, ncalls_mapped, ncalls_ref_allele_too_long;
  uint64_t nlines, nlines_too_long, nlines_match_ref, nlines_mapped;
  uint64_t nvars, nallele_too_long, nvars_printed; // decomposed ALTs
} DecomposeStats;

CallDecomp* call_decomp_init(FILE *fout);
void call_decomp_destroy(CallDecomp *dc);
scoring_t* call_decomp_get_scoring(CallDecomp *dc);

void call_decomp_cpy_stats(DecomposeStats *stats, const CallDecomp *dc);

void acall_decompose(CallDecomp *dc, const AlignedCall *call,
                     size_t max_line_len, size_t max_allele_len);

#endif /* ALIGNED_CALL_ */
