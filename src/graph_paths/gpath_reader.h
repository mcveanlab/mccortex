#ifndef GPATH_READER_H_
#define GPATH_READER_H_

#include "file_filter.h"
#include "db_graph.h"
#include "cJSON/cJSON.h"

typedef struct
{
  StrBuf hdrstr;
  cJSON *json;
  FileFilter fltr; // colour filter
  // Header is decomposed here
  // size_t kmer_size, num_paths, path_bytes, kmers_with_paths;
  size_t ncolours;
  cJSON **colours_json;
} GPathReader;

#include "objbuf_macro.h"
create_objbuf(gpfile_buf, GPathFileBuffer, GPathReader);

// 4MB buffer
// #define CTP_BUF_SIZE (4*ONE_MEGABYTE)

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GPathReader and returns 1
int gpath_reader_open(GPathReader *file, char *path, bool fatal);

// mode is "r", "r+" etc.
int gpath_reader_open2(GPathReader *file, char *path, const char *mode, bool fatal);

void gpath_reader_check(const GPathReader *file, size_t kmer_size, size_t ncols);
void gpath_reader_load(GPathReader *file, bool dont_add_kmers, dBGraph *db_graph);
void gpath_reader_close(GPathReader *file);

// Fetch information from header
size_t gpath_reader_get_kmer_size(const GPathReader *file);
size_t gpath_reader_get_num_kmers(const GPathReader *file);
size_t gpath_reader_get_num_paths(const GPathReader *file);
size_t gpath_reader_get_path_bytes(const GPathReader *file);
const char* gpath_reader_get_sample_name(const GPathReader *file, size_t idx);

// Copy sample names into the graph
void gpath_reader_load_sample_names(const GPathReader *file, dBGraph *db_graph);

// Get max mem required to load ctp files
void gpath_reader_max_mem_req(GPathReader *files, size_t nfiles,
                              size_t ncols, size_t graph_capacity,
                              bool store_nseen_klen,
                              bool split_lists, bool use_hash,
                              size_t *min_mem_ptr, size_t *max_mem_ptr);

#endif /* GPATH_READER_H_ */
