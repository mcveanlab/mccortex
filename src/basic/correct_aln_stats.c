#include "global.h"
#include "correct_aln_stats.h"

#define INIT_BUFLEN 1024

static inline void merge_arrays(uint64_t *restrict a,
                                const uint64_t *restrict b,
                                size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) a[i] += b[i];
}

void correct_aln_stats_alloc(CorrectAlnStats *stats)
{
  // Gap Histogram
  stats->histgrm_len = INIT_BUFLEN;
  stats->gap_ins_histgrm = calloc2(stats->histgrm_len, sizeof(uint64_t));
  stats->gap_err_histgrm = calloc2(stats->histgrm_len, sizeof(uint64_t));
}

void correct_aln_stats_dealloc(CorrectAlnStats *stats)
{
  free(stats->gap_ins_histgrm);
  free(stats->gap_err_histgrm);
}

void correct_aln_stats_zero(CorrectAlnStats *stats)
{
  memset(stats->gap_err_histgrm, 0, stats->histgrm_len * sizeof(uint64_t));
  memset(stats->gap_ins_histgrm, 0, stats->histgrm_len * sizeof(uint64_t));
  stats->num_gap_attempts = 0;
  stats->num_gap_successes = 0;
  stats->num_gaps_disagreed = 0;
  stats->num_gaps_too_short = 0;
}

void correct_aln_stats_merge(CorrectAlnStats *restrict dst,
                             CorrectAlnStats *restrict src)
{
  correct_aln_stats_cap(dst, src->histgrm_len);
  merge_arrays(dst->gap_err_histgrm, src->gap_err_histgrm, src->histgrm_len);
  merge_arrays(dst->gap_ins_histgrm, src->gap_ins_histgrm, src->histgrm_len);
  dst->num_gap_attempts += src->num_gap_attempts;
  dst->num_gap_successes += src->num_gap_successes;
  dst->num_gaps_disagreed += src->num_gaps_disagreed;
  dst->num_gaps_too_short += src->num_gaps_too_short;
}

void correct_aln_stats_cap(CorrectAlnStats *stats, size_t max_gap)
{
  if(stats->histgrm_len < max_gap) {
    max_gap = roundup2pow(max_gap);
    stats->gap_ins_histgrm = realloc2(stats->gap_ins_histgrm, max_gap*sizeof(uint64_t));
    stats->gap_err_histgrm = realloc2(stats->gap_err_histgrm, max_gap*sizeof(uint64_t));
    // Zero new memory
    size_t newmem = (max_gap-stats->histgrm_len)*sizeof(uint64_t);
    memset(stats->gap_ins_histgrm+stats->histgrm_len, 0, newmem);
    memset(stats->gap_err_histgrm+stats->histgrm_len, 0, newmem);
    stats->histgrm_len = max_gap;
  }
}
