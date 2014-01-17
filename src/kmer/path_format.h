#ifndef PATH_FORMAT_H_
#define PATH_FORMAT_H_

#include "graph_typedef.h"
#include "path_file_filter.h"
#include "path_store.h"

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
                           boolean fatal, const char *path);

// Get min number of colours needed to load the files
size_t paths_get_min_usedcols(PathFileReader *files, size_t num_files);

boolean paths_merge_needs_tmp(PathFileReader *files, size_t num_files);

// If tmppaths != NULL, do merge
// if insert is true, insert missing kmers into the graph
void paths_format_load(PathFileReader *file, dBGraph *db_graph,
                        boolean insert_missing_kmers);

// db_graph.pdata must be big enough to hold all this data or we exit
void paths_format_merge(PathFileReader *files, size_t num_files,
                         boolean insert_missing_kmers,
                         uint8_t *tmpdata, size_t tmpdatasize, dBGraph *db_graph);

//
// Write
//

// returns number of bytes written
size_t paths_format_write_header_core(const PathFileHeader *header, FILE *fout);
// returns number of bytes written
size_t paths_format_write_header(const PathFileHeader *header, FILE *fout);

void paths_format_write_optimised_paths(dBGraph *db_graph, FILE *fout);

#endif /* PATH_FORMAT_H_ */
