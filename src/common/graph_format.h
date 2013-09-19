#ifndef GRAPH_FORMAT_H_
#define GRAPH_FORMAT_H_

#include <inttypes.h>

#include "file_reader.h"
#include "db_graph.h"

// graph file format version
#define CTX_GRAPH_FILEFORMAT 6
// 4MB buffer
#define CTX_BUF_SIZE (4UL<<20)

typedef struct
{
  uint32_t version, kmer_size, num_of_bitfields, num_of_cols, max_col;
  uint64_t num_of_kmers;
  // Cleaning info etc for each colour
  GraphInfo *ginfo;
  uint32_t capacity; // number of colours malloc'd
} GraphFileHeader;

// Get an array of colour indices for a binary.
// arr should be of length num_of_cols
// in.c2.ctx {0,1}
// in.c2.ctx:1 {1}
uint32_t graph_file_get_ncols(const char *path, uint32_t max_col);
void graph_file_parse_colours(const char *str, uint32_t *arr, uint32_t max_col);

void graph_header_alloc(GraphFileHeader *header, size_t num_of_cols);
void graph_header_dealloc(GraphFileHeader *header);
void graph_header_cpy(GraphFileHeader *dst, const GraphFileHeader *src);

// If fatal == false, returns -1 on error
// path is used when reporting errors
int graph_file_read_header(FILE *fh, GraphFileHeader *header,
                           boolean fatal, const char *path);

// Returns number of bytes read
size_t graph_file_read_kmer(FILE *fh, GraphFileHeader *header, const char *path,
                            uint64_t *bkmer, Covg *covgs, Edges *edges);

// returns 0 if cannot read, 1 otherwise
boolean graph_file_probe(const char* ctx_path, boolean *valid_ctx,
                         GraphFileHeader *gheader);

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

uint32_t graph_load_colour(const char *path, const SeqLoadingPrefs *prefs,
                           SeqLoadingStats *stats, uint32_t colour);

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

//
// Writing
//

void graph_write_empty(const dBGraph *db_graph, FILE *fh, uint32_t num_of_cols);

// Dump a single colour into an existing binary
// FILE *fh must already point to the first bkmer
// if merge is true, read existing covg and edges and combine with outgoing
void graph_file_write_colour(dBGraph *db_graph, Colour graphcol,
                             Colour intocol, uint32_t num_of_cols, FILE *fh);

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
