#ifndef GPATH_READER_H_
#define GPATH_READER_H_

#include "file_filter.h"
#include "db_graph.h"
#include "cJSON/cJSON.h"

#include "common_buffers.h"

#define CTP_FORMAT_VERSION 4

typedef struct
{
  StreamBuffer strmbuf; 
  gzFile gz;

  // For parsing input
  StrBuf line;
  SizeBuffer numbuf;

  FileFilter fltr; // colour filter
  StrBuf hdrstr;
  cJSON *json;

  // Header is decomposed here
  // size_t kmer_size, num_paths, path_bytes, kmers_with_paths;
  int version;
  size_t ncolours;
  cJSON **colours_json;
} GPathReader;

#define GPATH_ADD_MISSING_KMERS   0
#define GPATH_DIE_MISSING_KMERS   1
#define GPATH_SKIP_MISSING_KMERS  2

#include "madcrowlib/madcrow_buffer.h"

madcrow_buffer(gpfile_buf, GPathFileBuffer, GPathReader);

/**
 Parse line with format:
  [FR] [njuncs] [nseen0,nseen1,...] [juncs:ACAGT] ([seq=] [juncpos=])?
 */
void link_line_parse(const StrBuf *line, int version, const FileFilter *fltr,
                     bool *fw, size_t *njuncs,
                     SizeBuffer *counts, StrBuf *juncs,
                     StrBuf *seq, SizeBuffer *juncpos);

// Reads line <kmer> <num_links>
// Calls die() on error
// Returns true unless end of file
bool gpath_reader_read_kmer(GPathReader *file, StrBuf *kmer, size_t *num_links);

// Reads line [FR] <num_links>
// Calls die() on error
// Returns true unless end of link entries
bool gpath_reader_read_link(GPathReader *file,
                            bool *fw, size_t *njuncs,
                            SizeBuffer *countbuf, StrBuf *juncs,
                            StrBuf *seq, SizeBuffer *juncpos);


// Open file, exits on error
// if successful creates a new GPathReader
void gpath_reader_open(GPathReader *file, const char *path);

// mode is "r", "r+" etc.
// If there is no 'into' filter (e.g. 0,1:in.ctx 0,1 is 'into filter'),
// load into `into_offset..into_offset+N-1`
void gpath_reader_open2(GPathReader *file, const char *path, const char *mode,
                        size_t into_offset);

void gpath_reader_check(const GPathReader *file, size_t kmer_size, size_t ncols);

// @kmer_flags must be one of:
//   GPATH_ADD_MISSING_KMERS - add kmers to the graph before loading path
//   GPATH_DIE_MISSING_KMERS - die with error if cannot find kmer
//   GPATH_SKIP_MISSING_KMERS - skip paths where kmer is not in graph
void gpath_reader_load(GPathReader *file, int kmer_flags, dBGraph *db_graph);
void gpath_reader_close(GPathReader *file);

//
// Reading without loading into a graph
//

// Reads line <kmer> <num_links>
// Calls die() on error
// Returns true unless end of file
bool gpath_reader_read_kmer(GPathReader *file, StrBuf *kmer, size_t *num_links);

// Reads line [FR] <num_links>
// Calls die() on error
// Returns true unless end of link entries
bool gpath_reader_read_link(GPathReader *file,
                            bool *fw, size_t *njuncs,
                            SizeBuffer *countbuf, StrBuf *juncs,
                            StrBuf *seq, SizeBuffer *juncpos);


//
// Fetch information from header
//
size_t gpath_reader_get_kmer_size(const GPathReader *file);
size_t gpath_reader_get_num_kmers(const GPathReader *file);
size_t gpath_reader_get_num_paths(const GPathReader *file);
size_t gpath_reader_get_path_bytes(const GPathReader *file);
const char* gpath_reader_get_sample_name(const GPathReader *file, size_t idx);

// Copy sample names into the graph
void gpath_reader_load_sample_names(const GPathReader *file, dBGraph *db_graph);

void gpath_reader_load_contig_hist(cJSON *json_root, const char *path,
                                   size_t fromcol, ZeroSizeBuffer *hist);

//
// Memory Calculations
//

// Get max mem required to load ctp files
void gpath_reader_max_mem_req(GPathReader *files, size_t nfiles,
                              size_t ncols, size_t graph_capacity,
                              bool store_nseen,
                              bool split_lists, bool use_hash,
                              size_t *min_mem_ptr, size_t *max_mem_ptr);

// Returns sum of memory
// sets @max_file_mem_ptr to the max for a single file memory
size_t gpath_reader_sum_mem(GPathReader *files, size_t nfiles,
                            size_t ncols, bool count_nseen, bool use_gphash,
                            size_t *max_file_mem_ptr);

size_t gpath_reader_mem_req(GPathReader *files, size_t nfiles,
                            size_t ncols, size_t max_mem,
                            bool count_nseen);

void gpath_reader_alloc_gpstore(GPathReader *files, size_t nfiles,
                                size_t mem, bool count_nseen,
                                dBGraph *db_graph);

#endif /* GPATH_READER_H_ */
