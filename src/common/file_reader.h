#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <string_buffer.h>
#include "seq_file.h"

#include "db_graph.h"

// Used in parse_filelist
#define READ_FALIST 0
#define READ_COLOURLIST 1

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
  unsigned long kmers_loaded, unique_kmers;
  // Used for binaries and colourlists
  unsigned long num_of_colours_loaded;
} SeqLoadingStats;

// Functions for dealing with file loading statistics
SeqLoadingStats* seq_loading_stats_create(unsigned long readlen_arrsize);
void seq_loading_stats_sum(SeqLoadingStats* dst, SeqLoadingStats* src);
void seq_loading_stats_free(SeqLoadingStats* stats);

// Stucture for specifying how to load data
typedef struct
{
  Colour into_colour;
  boolean merge_colours; // Load all data into only one colour
  boolean boolean_covgs; // Update covg by at most 1

  // loading sequence
  boolean load_seq;
  char quality_cutoff, ascii_fq_offset;
  int homopolymer_cutoff;
  boolean remove_dups_se, remove_dups_pe;

  // loading binaries
  boolean load_binaries;
  boolean must_exist_in_graph;
  // if empty_colours is true an error is thrown if a kmer is repeated
  boolean empty_colours;

  boolean update_ginfo;
  dBGraph *db_graph;
} SeqLoadingPrefs;


void parse_filelists(const char *list_path1, const char *list_path2,
                     uint8_t are_colour_lists,
                     SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                     void (*read_func)(read_t *r1, read_t *r2,
                                       int qoffset1, int qoffset2,
                                       SeqLoadingPrefs *prefs,
                                       SeqLoadingStats *stats, void *ptr),
                     void *reader_ptr);

// Given a 'filelist' file, check all files pointed to exist.
// If path_array is not NULL, populate it with the paths.
// Exits with error if a file doesn't exist or isn't readable.
// Return the number of files pointed to.
// ginfo is optional and is used to update sample names
uint32_t load_paths_from_filelist(const char *filelist_path, char **path_array,
                                  boolean sample_names_permitted,
                                  GraphInfo *ginfo, boolean *has_sample_names);

// filename is a list of files, one for each colour (with optional second column
// of sample-ids). Check they all exists, there are not too many, and that each
// of them contains a list of valid binaries.

// Returns number for files in file passed
// Check that ctxlists contain only valid ctx files
// Check that colourlists don't contain any ctx files
// `colours` is only updated for colour lists
uint32_t check_colour_or_ctx_list(const char *list_path, boolean is_colourlist,
                                  boolean binaries_allowed, boolean seq_allowed,
                                  uint32_t kmer_size, uint32_t *colours);

#endif /* FILE_READER_H_ */
