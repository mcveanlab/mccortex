#ifndef LOADING_STATS_H_
#define LOADING_STATS_H_

// Stucture for statistics on loading sequence and cortex binary files
typedef struct
{
  size_t se_colourlists_loaded, pe_colourlist_pairs_loaded;
  size_t se_filelists_loaded, pe_filelist_pairs_loaded;
  size_t num_files_loaded, binaries_loaded;

  size_t num_se_reads, num_pe_reads;
  size_t total_good_reads, total_bad_reads, total_dup_reads;
  size_t total_bases_read, total_bases_loaded;
  size_t *readlen_count_array;
  size_t readlen_count_array_size;
  size_t contigs_loaded, kmers_loaded, unique_kmers;
  // Used for binaries and colourlists
  size_t num_of_colours_loaded;
} SeqLoadingStats;

// Functions for dealing with file loading statistics
SeqLoadingStats* seq_loading_stats_create(size_t readlen_arrsize);
void seq_loading_stats_sum(SeqLoadingStats* dst, SeqLoadingStats* src);
void seq_loading_stats_free(SeqLoadingStats* stats);

#endif /* LOADING_STATS_H_ */
