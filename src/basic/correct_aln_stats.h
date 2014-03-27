#ifndef CORRECT_ALN_STATS_H_
#define CORRECT_ALN_STATS_H_

typedef struct
{
  uint64_t *gap_ins_histgrm, *gap_err_histgrm, histgrm_len;
  // Count cases where traversal worked but paths disagreed with remaining
  // contig
  uint64_t num_gap_attempts, num_gap_successes;
  uint64_t num_gaps_disagreed, num_gaps_too_short;
} CorrectAlnStats;

void correct_aln_stats_alloc(CorrectAlnStats *stats);
void correct_aln_stats_dealloc(CorrectAlnStats *stats);
void correct_aln_stats_zero(CorrectAlnStats *stats);
void correct_aln_stats_merge(CorrectAlnStats *restrict dst,
                             CorrectAlnStats *restrict src);
void correct_aln_stats_cap(CorrectAlnStats *stats, size_t max_gap);

// Save gap size distribution
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void correct_aln_stats_dump(const char *path,
                            const uint64_t *arr, size_t arrlen,
                            size_t kmer_size, bool insert_sizes,
                            size_t nreads);

// DEV: add functions for recording ins gap, err gap

#endif /* CORRECT_ALN_STATS_H_ */
