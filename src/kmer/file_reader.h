#ifndef FILE_READER_H_
#define FILE_READER_H_

#include "seq_file.h"
#include "graph_typedef.h"

// Stucture for specifying how to load data
typedef struct
{
  dBGraph *db_graph;

  // Graphs
  // boolean merge_colours; // Load all data into only one colour
  boolean boolean_covgs; // Update covg by at most 1
  boolean must_exist_in_graph;
  // if empty_colours is true an error is thrown a kmer from a binary is
  // already in the graph
  boolean empty_colours;

  // Sequence
  Colour into_colour;
  char quality_cutoff, ascii_fq_offset;
  int homopolymer_cutoff;
  boolean remove_dups_se, remove_dups_pe;

} SeqLoadingPrefs;

// Stucture for statistics on loading sequence and cortex binary files
typedef struct
{
  unsigned long se_colourlists_loaded, pe_colourlist_pairs_loaded;
  unsigned long se_filelists_loaded, pe_filelist_pairs_loaded;
  unsigned long num_files_loaded, binaries_loaded;

  unsigned long num_se_reads, num_pe_reads;
  unsigned long total_good_reads, total_bad_reads, total_dup_reads;
  unsigned long total_bases_read, total_bases_loaded;
  unsigned long *readlen_count_array;
  unsigned long readlen_count_array_size;
  unsigned long contigs_loaded, kmers_loaded, unique_kmers;
  // Used for binaries and colourlists
  unsigned long num_of_colours_loaded;
} SeqLoadingStats;

// Functions for dealing with file loading statistics
SeqLoadingStats* seq_loading_stats_create(size_t readlen_arrsize);
void seq_loading_stats_sum(SeqLoadingStats* dst, SeqLoadingStats* src);
void seq_loading_stats_free(SeqLoadingStats* stats);

#endif /* FILE_READER_H_ */
