#ifndef PATH_FORMAT_H_
#define PATH_FORMAT_H_

#include "graph_typedef.h"

// file format version
#define CTX_PATH_FILEFORMAT 1

void paths_format_write(const dBGraph *db_graph, const binary_paths_t *paths,
                        const char *path);

// if insert is true, insert missing kmers into the graph
void paths_format_read(dBGraph *db_graph, binary_paths_t *paths,
                       boolean insert, const char *path);

// Returns false if cannot read otherwise true
boolean paths_format_probe(const char *path, boolean *valid_paths_file,
                           uint32_t *kmer_size_ptr, uint32_t *num_of_cols_ptr,
                           uint64_t *num_paths_ptr, uint64_t *num_path_bytes_ptr,
                           uint64_t *num_path_kmers_ptr);

void paths_format_filename(const char *path, char *out);

#endif
