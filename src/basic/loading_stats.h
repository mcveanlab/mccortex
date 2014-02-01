#ifndef LOADING_STATS_H_
#define LOADING_STATS_H_

// Stucture for statistics on loading sequence and cortex binary files
typedef struct
{
  size_t num_files_loaded, ctx_files_loaded;
  size_t num_se_reads, num_pe_reads;
  size_t total_good_reads, total_bad_reads, total_dup_reads;
  size_t total_bases_read, total_bases_loaded;
  size_t contigs_loaded, kmers_loaded, unique_kmers;
  size_t num_of_colours_loaded; // ctx files only
} LoadingStats;

#define LOAD_STATS_INIT_MACRO { \
  .num_files_loaded = 0, .ctx_files_loaded = 0, \
  .num_se_reads = 0, .num_pe_reads = 0, \
  .total_good_reads = 0, .total_bad_reads = 0, .total_dup_reads = 0, \
  .total_bases_read = 0, .total_bases_loaded = 0, \
  .contigs_loaded = 0, .kmers_loaded = 0, .unique_kmers = 0, \
  .num_of_colours_loaded = 0 \
}

// Functions for dealing with file loading statistics
void loading_stats_init(LoadingStats *stats);
void loading_stats_merge(LoadingStats* dst, const LoadingStats* src);

#endif /* LOADING_STATS_H_ */
