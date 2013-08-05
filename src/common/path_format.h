#ifndef PATH_FORMAT_H_
#define PATH_FORMAT_H_

#include "graph_typedef.h"

// path file format version
#define CTX_PATH_FILEFORMAT 1

typedef struct
{
  uint32_t version, kmer_size, num_of_cols;
  uint64_t num_of_paths, num_path_bytes, num_kmers_with_paths;
  StrBuf *sample_names;
  size_t capacity; // how many sample names have been malloc'd
} PathFileHeader;

// Path File Format:
// -- Header --
// "PATHS"<uint32_t:version><uint32_t:kmersize><uint32_t:num_of_cols>
// <uint64_t:num_of_paths><uint64_t:num_path_bytes><uint64_t:num_kmers_with_paths>
// -- Colours --
// <uint32_t:sname_len><uint8_t x sname_len:sample_name> x num_of_cols
// <uint64_t:num_inferred_kmers> x num_of_cols
// -- Inferred kmers --
// [<binarykmer> x num_inferred_kmers[col]] x num_of_cols
// -- Data --
// <uint8_t:path_data>
// <binarykmer><uint64_t:path_index_fw><uint64_t:path_index_rv>

void paths_format_filename(const char *file_path, char *out);

void paths_header_alloc(PathFileHeader *header, size_t num_of_cols);
void paths_header_dealloc(PathFileHeader *header);

// Returns number of bytes read or -1 on error (if fatal == false)
int paths_file_read_header(FILE *fh, PathFileHeader *header,
                           boolean fatal, const char *path);

// Returns false if cannot read otherwise true
boolean paths_file_probe(const char *ctp_path, boolean *valid_ctp,
                         PathFileHeader *pheader);

// Check a header and graph are compatible
void paths_graph_compatibility_check(const PathFileHeader *pheader,
                                     const dBGraph *db_graph);

// If tmppaths != NULL, do merge
// if insert is true, insert missing kmers into the graph
void paths_format_merge(const char *path, PathFileHeader *pheader,
                        dBGraph *db_graph, PathStore *paths,
                        PathStore *tmppaths, boolean insert_missing_kmers);

void paths_format_read(const char *path, PathFileHeader *pheader,
                       dBGraph *db_graph, PathStore *paths,
                       boolean insert_missing_kmers);

void paths_format_write(const dBGraph *db_graph, const PathStore *paths,
                        PathFileHeader *header, const char *file_path);

#endif
