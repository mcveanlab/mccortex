#ifndef GPATH_READER_H_
#define GPATH_READER_H_

#include "file_filter.h"
#include "db_graph.h"
#include "cJSON/cJSON.h"

typedef struct
{
  gzFile gz;
  FileFilter fltr; // colour filter
  StrBuf hdrstr;
  cJSON *json;
  // Header is decomposed here
  // size_t kmer_size, num_paths, path_bytes, kmers_with_paths;
  size_t ncolours;
  cJSON **colours_json;
} GPathReader;

#define GPATH_ADD_MISSING_KMERS   0
#define GPATH_DIE_MISSING_KMERS   1
#define GPATH_SKIP_MISSING_KMERS  2

#include "objbuf_macro.h"
create_objbuf(gpfile_buf, GPathFileBuffer, GPathReader);

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

// Fetch information from header
size_t gpath_reader_get_kmer_size(const GPathReader *file);
size_t gpath_reader_get_num_kmers(const GPathReader *file);
size_t gpath_reader_get_num_paths(const GPathReader *file);
size_t gpath_reader_get_path_bytes(const GPathReader *file);
const char* gpath_reader_get_sample_name(const GPathReader *file, size_t idx);

// Copy sample names into the graph
void gpath_reader_load_sample_names(const GPathReader *file, dBGraph *db_graph);

//
// Memory Calculations
//

// Get max mem required to load ctp files
void gpath_reader_max_mem_req(GPathReader *files, size_t nfiles,
                              size_t ncols, size_t graph_capacity,
                              bool store_nseen_klen,
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

// Create a path store that does not tracks path counts
void gpath_reader_alloc_gpstore(GPathReader *files, size_t nfiles,
                                size_t mem, bool count_nseen,
                                dBGraph *db_graph);

#endif /* GPATH_READER_H_ */
