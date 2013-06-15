#ifndef GRAPH_TYPEDEF_H_
#define GRAPH_TYPEDEF_H_

#include "string_buffer.h"
#include "hash_table.h"

typedef uint8_t Colour;
typedef uint8_t Edges;
typedef uint32_t Covg;

#define COVG_MAX UINT_MAX

typedef enum
{
  forward = 0,
  reverse = 1
} Orientation;

typedef BinaryKmerPtr Key;

typedef struct {
  uint8_t *const store, *const end;
  const size_t size;
  uint8_t *next;
  size_t num_paths;
} binary_paths_t;

typedef struct
{
  boolean tip_clipping, remv_low_cov_sups, remv_low_cov_nodes;

  // Thesholds are zero if not used (e.g. remv_low_cov_sups == false)
  Covg remv_low_cov_sups_thresh, remv_low_cov_nodes_thresh;

  // Cleaning a low covg sample against cleaned pool of population
  boolean cleaned_against_another_graph;
  StrBuf* cleaned_against_graph_name;
} ErrorCleaning;

typedef struct
{
  StrBuf *sample_names[NUM_OF_COLOURS];
  uint64_t total_sequence[NUM_OF_COLOURS];
  uint32_t mean_read_length[NUM_OF_COLOURS];
  long double seq_err[NUM_OF_COLOURS];
  ErrorCleaning cleaning[NUM_OF_COLOURS];
  uint32_t num_of_colours_loaded, num_of_shades_loaded;
} GraphInfo;

typedef struct
{
  HashTable ht;
  uint32_t kmer_size;
  // Optional fields:
  Covg (*covgs)[NUM_OF_COLOURS]; // [hkey][col]
  Edges *edges;
  Edges (*col_edges)[NUM_OF_COLOURS]; // [hkey][col]
  uint64_t *node_in_cols[NUM_OF_COLOURS]; // [col][hkey/64]
  // path data
  uint64_t *kmer_paths;
  binary_paths_t pdata;
  // Loading reads
  uint64_t *readstrt;
  // Info stored here:
  GraphInfo ginfo;
} dBGraph;

#endif
