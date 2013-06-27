#ifndef BINARY_FORMAT_H_
#define BINARY_FORMAT_H_

#include <inttypes.h>

#include "file_reader.h"
#include "db_graph.h"

extern const int CURR_CTX_VERSION;
extern const char CTX_MAGIC_WORD[7];

typedef struct
{
  uint32_t version, kmer_size, num_of_bitfields, num_of_colours;
  uint64_t num_of_kmers, *total_seq_loaded;
  uint32_t *mean_read_lengths;
  StrBuf **sample_names;
  long double *seq_err_rates;
  ErrorCleaning **err_cleaning;
} BinaryFileHeader;


size_t binary_read_header(FILE *fh, BinaryFileHeader *header, const char *path);
size_t binary_read_kmer(FILE *fh, BinaryFileHeader *header, const char *path,
                        uint64_t *bkmer, Covg *covgs, Edges *edges);

// After calling binary_read_header you must call:
void binary_header_destroy(BinaryFileHeader *header);

void binary_write_header(FILE *fh, const BinaryFileHeader *header);
void binary_write_kmer(FILE *fh, const BinaryFileHeader *h,
                       const uint64_t *bkmer, const Covg *covgs,
                       const Edges *edges);

// returns 0 if cannot read, 1 otherwise
char binary_probe(const char* path, boolean *is_ctx,
                  uint32_t *kmer_size_ptr, uint32_t *num_of_colours_ptr,
                  uint64_t *num_of_kmers);

// Load a binary, putting the first colour into `load_first_colour_into`
//
// if only_load_if_in_colour is >= 0 only kmers with coverage in existing
//    colour only_load_if_in_colour will be loaded.
// if empty_colours != 0 an error is thrown if a node already exists
// if load_as_union != 0 then we only increment covg if it is zero
uint32_t binary_load(const char *path, dBGraph *graph,
                     SeqLoadingPrefs *prefs, SeqLoadingStats *stats);

// This function will dump valid binaries by not printing edges to nodes that
// are not themselves printed
// If you don't want to/care about graph_info, pass in NULL
// If you want to print all nodes pass condition as NULL
// start_col is ignored unless colours is NULL
uint64_t binary_dump_graph(const char *path, dBGraph *graph,
                           uint32_t version,
                           const Colour *colours, Colour start_col,
                           uint32_t num_of_cols);

#endif /* BINARY_FORMAT_H_ */
