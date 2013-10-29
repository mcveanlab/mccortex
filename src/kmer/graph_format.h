#ifndef GRAPH_FORMAT_H_
#define GRAPH_FORMAT_H_

#include <inttypes.h>

#include "file_reader.h"
#include "db_graph.h"
#include "graph_file_filter.h"

// graph file format version
#define CTX_GRAPH_FILEFORMAT 6

// Get an array of colour indices for a binary.
// arr should be of length num_of_cols
// in.c2.ctx {0,1}
// in.c2.ctx:1 {1}
uint32_t graph_file_get_ncols(const char *path, uint32_t max_col);
void graph_file_parse_colours(const char *str, uint32_t *arr, uint32_t max_col);

void graph_header_alloc(GraphFileHeader *header, size_t num_of_cols);
void graph_header_dealloc(GraphFileHeader *header);

// Copy non-colour specific values
void graph_header_global_cpy(GraphFileHeader *dst, const GraphFileHeader *src);

// If fatal == false, returns -1 on error
// path is used when reporting errors
int graph_file_read_header(FILE *fh, GraphFileHeader *header,
                           boolean fatal, const char *path);

// Returns number of bytes read
size_t graph_file_read_kmer(FILE *fh, const GraphFileHeader *h, const char *path,
                            uint64_t *bkmer, Covg *covgs, Edges *edges);

// GFF: unproposed
// size_t graph_file_read_kmer2(const GraphFileFilter *gff,
//                              uint64_t *bkmer, Covg *covgs, Edges *edges);

// returns 0 if cannot read, 1 otherwise
boolean graph_file_probe(const char* ctx_path, boolean *valid_ctx,
                         GraphFileHeader *gheader);

// NOT PROPOSED:
// boolean graph_file_probe2(const char* ctx_path, boolean *valid_ctx,
//                           GraphFileHeader *gheader, GraphFileFilter *gff);

// if only_load_if_in_colour is >= 0, only kmers with coverage in existing
// colour only_load_if_in_colour will be loaded.
// if clean_colours != 0 an error is thrown if a node already exists
// returns the number of colours in the binary
// If stats != NULL, updates:
//   stats->num_of_colours_loaded
//   stats->kmers_loaded
//   stats->total_bases_read
//   stats->binaries_loaded
// If header is != NULL, header will be stored there.  Be sure to free.
uint32_t graph_load(const char *path,
                    const SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                    GraphFileHeader *header);

// PROPOSED:
// closes file but does not dealloc
void graph_load2(GraphFileReader *file, const SeqLoadingPrefs *prefs,
                 SeqLoadingStats *stats);

// propose delete
uint32_t graph_load_colour(const char *path, const SeqLoadingPrefs *prefs,
                           SeqLoadingStats *stats, uint32_t colour);

// Load a kmer and write to a file one kmer at a time
// Optionally filter a against the graph currently loaded
//   (i.e. only keep nodes and edges that are in the graph)
// Same functionality as graph_files_merge, but faster if dealing with only one
// input file. Reads in and dumps one kmer at a time
// parameter: flatten: if true merge colours into one
size_t graph_stream_filter(const char *out_ctx_path, char *in_ctx_path,
                           boolean flatten, const dBGraph *db_graph,
                           const char *intersect_gname);

// PROPOSED: written
size_t graph_stream_filter2(const char *out_ctx_path, GraphFileReader *file,
                            const dBGraph *db_graph, const char *intersect_gname);

// ctx_num_cols and ctx_max_cols are the numbers returned from binary_probe
// if merge pool colour 0 from each binary into colour 0, 1 -> 1 etc.
// if flatten, pool all colours into colour 0
// if intersect only load kmers that are already in the hash table
// returns number of kmers written
size_t graph_files_merge(char *out_ctx_path, char **binary_paths,
                         size_t num_binaries,
                         uint32_t ctx_num_cols[num_binaries],
                         uint32_t ctx_max_cols[num_binaries],
                         boolean merge, boolean flatten,
                         const char *intersect_gname,
                         dBGraph *db_graph);

// PROPOSED:
size_t graph_files_merge2(char *out_ctx_path, GraphFileReader *files,
                          size_t num_files,
                          const char *intersect_gname, dBGraph *db_graph);

//
// Writing
//

void graph_write_empty(const dBGraph *db_graph, FILE *fh, size_t num_of_cols);

// Dump a single colour into an existing binary
// FILE *fh must already point to the first bkmer
// if merge is true, read existing covg and edges and combine with outgoing
void graph_file_write_colour(const dBGraph *db_graph, Colour graphcol,
                             Colour intocol, size_t file_ncols, FILE *fh);

void graph_file_write_colours(const dBGraph *db_graph, Colour graphcol,
                              Colour intocol, size_t write_ncols,
                              size_t file_ncols, FILE *fh);

// Returns number of bytes written
size_t graph_write_header(FILE *fh, const GraphFileHeader *header);
size_t graph_write_kmer(FILE *fh, const GraphFileHeader *h,
                        const uint64_t *bkmer, const Covg *covgs,
                        const Edges *edges);

// If you don't want to/care about graph_info, pass in NULL
// If you want to print all nodes pass condition as NULL
// start_col is ignored unless colours is NULL
// returns number of nodes dumped
uint64_t graph_file_save(const char *path, dBGraph *graph,
                         uint32_t version,
                         const Colour *colours, Colour start_col,
                         uint32_t num_of_cols);

void graph_write_status(uint64_t nkmers, size_t ncols,
                        const char *path, uint32_t version);

#endif /* GRAPH_FORMAT_H_ */
