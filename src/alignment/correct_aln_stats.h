#ifndef CORRECT_ALN_STATS_H_
#define CORRECT_ALN_STATS_H_

#include "loading_stats.h"

#define ALN_STATS_MAX_GAP 128
#define ALN_STATS_MAX_FRAGLEN 1024

typedef struct
{
  // fraglen_histgrm is in BASES not kmers
  size_t fraglen_histgrm[ALN_STATS_MAX_FRAGLEN];
  // [exp_gap][act_gap]
  size_t gap_err_histgrm[ALN_STATS_MAX_GAP][ALN_STATS_MAX_GAP];
  // Count cases where traversal worked but paths disagreed with remaining contig
  uint64_t num_gap_attempts, num_gap_successes;
  uint64_t num_paths_disagreed, num_gaps_too_short;
} CorrectAlnStats;

typedef struct {
  uint32_t gap_len;
  bool traversed, paths_disagreed, gap_too_short;
} TraversalResult;

void correct_aln_stats_reset(CorrectAlnStats *stats);
void correct_aln_stats_merge(CorrectAlnStats *restrict dst,
                             CorrectAlnStats *restrict src);

static inline void correct_aln_stats_update(CorrectAlnStats *stats,
                                            TraversalResult result)
{
  stats->num_gap_attempts++;
  stats->num_gap_successes += result.traversed;
  stats->num_paths_disagreed += result.paths_disagreed;
  stats->num_gaps_too_short += result.gap_too_short;
}

// Sequencing error gap
void correct_aln_stats_add(CorrectAlnStats *stats,
                           size_t exp_seq_gap, size_t act_gap);

// @exp_seq_gap does not include mate pair gap
void correct_aln_stats_add_mp(CorrectAlnStats *stats,
                              size_t exp_seq_gap, size_t gap_kmers,
                              size_t r1bases, size_t r2bases,
                              size_t kmer_size);

// Save gap size distribution
void correct_aln_stats_dump_gaps(const CorrectAlnStats *stats, const char *path);
// Save fragment size vector
void correct_aln_stats_dump_fraglen(const CorrectAlnStats *stats, const char *path);

void correct_aln_stats_print_summary(const CorrectAlnStats *stats,
                                     size_t num_reads, size_t num_read_pairs);

// Print summary stats and write output files
// @ht_num_kmers is the number of kmers loaded into the graph
void correct_aln_dump_stats(const LoadingStats *stats,
                            const CorrectAlnStats *gapstats,
                            const char *dump_seqgap_hist_path,
                            const char *dump_fraglen_hist_path,
                            size_t ht_num_kmers);

#endif /* CORRECT_ALN_STATS_H_ */
