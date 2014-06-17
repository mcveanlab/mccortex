#ifndef PATH_FORMAT_H_
#define PATH_FORMAT_H_

#include "cortex_types.h"
#include "path_file_reader.h"
#include "path_store.h"
#include "db_graph.h"
#include "db_node.h"

// path file format version
#define CTX_PATH_FILEFORMAT 1

// Path File Format:
// -- Header --
// "PATHS"<uint32_t:version><uint32_t:kmersize><uint32_t:num_of_cols>
// <uint64_t:num_of_paths><uint64_t:num_path_bytes><uint64_t:num_kmers_with_paths>
// -- Colours --
// <uint32_t:sname_len><uint8_t x sname_len:sample_name> x num_of_cols
// -- Data --
// <uint8_t:path_data>
// <binarykmer><uint64_t:path_index_fw><uint64_t:path_index_rv>

void paths_header_alloc(PathFileHeader *header, size_t num_of_cols);
void paths_header_dealloc(PathFileHeader *header);

// Set path header variables based on PathStore
void paths_header_update(PathFileHeader *header, const PathStore *paths);

// Returns number of bytes read or -1 on error (if fatal == false)
int paths_file_read_header(FILE *fh, PathFileHeader *header,
                           bool fatal, const char *path);

// Get min number of colours needed to load the files
size_t paths_get_max_usedcols(PathFileReader *files, size_t num_files);

// Get min memory required to load files. Returns memory required in bytes.
// remove_substr requires extra memory if only loading one file
// (as if loading two files)
size_t path_files_mem_required(const PathFileReader *files, size_t num_files,
                               bool remove_substr, bool use_path_hash,
                               size_t num_pcols, size_t extra_bytes);

// if insert is true, insert missing kmers into the graph
void paths_format_load(PathFileReader *file, bool insert_missing_kmers,
                       dBGraph *db_graph);

// Only load a given colour
// colour_idx is the index of an already specified colour
// "in.ctp:0,3,9,2" colour_idx=2 loads: "in.ctp:9"
void paths_load_colour(PathFileReader *pfile,
                       bool insert_missing_kmers,
                       size_t colour_idx, size_t intocol,
                       dBGraph *db_graph);

// Load 1 or more path files; can be called consecutively
// if `rmv_redundant` is true we remove non-informative paths
//  `thread_limit` is the number of threads to use for removing redundant paths
void paths_format_merge(PathFileReader *files, size_t num_files,
                        bool insert_missing_kmers,
                        bool rmv_redundant, size_t thread_limit,
                        dBGraph *db_graph);

//
// Write
//

// returns number of bytes written
size_t paths_format_write_header_core(const PathFileHeader *header, FILE *fout);
// returns number of bytes written
size_t paths_format_write_header(const PathFileHeader *header, FILE *fout);

void paths_format_write_optimised_paths_only(dBGraph *db_graph, FILE *fout);
void paths_format_write_optimised_paths(dBGraph *db_graph, FILE *fout);

#endif /* PATH_FORMAT_H_ */
