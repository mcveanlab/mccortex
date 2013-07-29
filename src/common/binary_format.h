#ifndef BINARY_FORMAT_H_
#define BINARY_FORMAT_H_

#include <inttypes.h>

#include "file_reader.h"
#include "db_graph.h"

extern const int CURR_CTX_VERSION;
extern const char CTX_MAGIC_WORD[7];

typedef struct
{
  uint32_t version, kmer_size, num_of_bitfields, num_of_cols;
  uint64_t num_of_kmers;
  // Cleaning info etc for each colour
  GraphInfo *ginfo;
  uint32_t capacity; // number of colours malloc'd
} BinaryFileHeader;

// Get an array of colour indices for a binary.
// arr should be of length num_of_cols
// in.c2.ctx {0,1}
// in.c2.ctx:1 {1}
uint32_t binary_get_num_colours(const char *path, uint32_t max_col);
void binary_parse_colour_array(const char *str, uint32_t *arr, uint32_t max_col);

// returns 0 if cannot read, 1 otherwise
char binary_probe(const char* path, boolean *is_ctx,
                  uint32_t *kmer_size_ptr, uint32_t *num_of_colours_ptr,
                  uint32_t *max_col_index, uint64_t *num_of_kmers);

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
uint32_t binary_load(const char *path, dBGraph *graph,
                     const SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                     BinaryFileHeader *header);

uint32_t binary_load_colour(const char *path, dBGraph *db_graph,
                            SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                            uint32_t colour);

size_t binary_read_header(FILE *fh, BinaryFileHeader *header, const char *path);

size_t binary_read_kmer(FILE *fh, BinaryFileHeader *header, const char *path,
                        uint64_t *bkmer, Covg *covgs, Edges *edges);

void binary_header_alloc(BinaryFileHeader *header, size_t num_of_cols);
void binary_header_realloc(BinaryFileHeader *header, size_t num_of_cols);
void binary_header_dealloc(BinaryFileHeader *header);
void binary_read_cpy_basic(BinaryFileHeader *dst, BinaryFileHeader *src);

void dump_empty_binary(dBGraph *db_graph, FILE *fh, uint32_t num_of_cols);

// Returns number of bytes written
size_t binary_write_header(FILE *fh, const BinaryFileHeader *header);
size_t binary_write_kmer(FILE *fh, const BinaryFileHeader *h,
                         const uint64_t *bkmer, const Covg *covgs,
                         const Edges *edges);

// If you don't want to/care about graph_info, pass in NULL
// If you want to print all nodes pass condition as NULL
// start_col is ignored unless colours is NULL
// returns number of nodes dumped
uint64_t binary_dump_graph(const char *path, dBGraph *graph,
                           uint32_t version,
                           const Colour *colours, Colour start_col,
                           uint32_t num_of_cols);

// Dump a single colour into an existing binary
// FILE *fh must already point to the first bkmer
// if merge is true, read existing covg and edges and combine with outgoing
void binary_dump_colour(dBGraph *db_graph, Colour graphcol,
                        Colour intocol, uint32_t num_of_cols, FILE *fh);

#endif /* BINARY_FORMAT_H_ */
