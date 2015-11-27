#ifndef SEQ_LOADING_STATS_H_
#define SEQ_LOADING_STATS_H_

// Stucture for statistics on loading sequence and cortex binary files
typedef struct
{
  // num_se_reads includes good reads, bad reads and duplicates etc.
  size_t num_se_reads, num_pe_reads;
  size_t num_good_reads, num_bad_reads, num_dup_se_reads, num_dup_pe_pairs;
  size_t total_bases_read, total_bases_loaded;
  size_t contigs_parsed, num_kmers_parsed, num_kmers_loaded, num_kmers_novel;
  uint64_t *col_nkmers, *col_sum_covgs;
  size_t ncols; // max number of colours loaded
} SeqLoadingStats;

#define SEQ_LOADING_STATS_INIT (SeqLoadingStats){ \
  .num_se_reads     = 0, .num_pe_reads       = 0, \
  .num_good_reads   = 0, .num_bad_reads      = 0, \
  .num_dup_se_reads = 0, .num_dup_pe_pairs   = 0, \
  .total_bases_read = 0, .total_bases_loaded = 0, \
  .contigs_parsed   = 0, .num_kmers_parsed   = 0, \
  .num_kmers_loaded = 0, .num_kmers_novel    = 0, \
  .col_nkmers = NULL, .col_sum_covgs = NULL, \
  .ncols = 0 \
}

// Functions for dealing with file loading statistics
#define seq_loading_stats_init(s) memset(s, 0, sizeof(SeqLoadingStats))
void seq_loading_stats_merge(SeqLoadingStats *dst, const SeqLoadingStats *src);

// @ht_num_kmers is the number of kmers loaded into the graph
void seq_loading_stats_print(const SeqLoadingStats *stats, size_t ht_num_kmers);

#endif /* SEQ_LOADING_STATS_H_ */
