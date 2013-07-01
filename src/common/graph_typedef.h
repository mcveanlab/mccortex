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
  const uint32_t num_of_cols;
  uint8_t *next;
  size_t num_paths;
} binary_paths_t;

// Thesholds are zero if not used (e.g. remv_low_cov_sups == false)
// cleaned_against_another_graph is for cleaning a low covg sample against
// cleaned pool of population
typedef struct
{
  boolean tip_clipping, remv_low_cov_sups, remv_low_cov_nodes;
  Covg remv_low_cov_sups_thresh, remv_low_cov_nodes_thresh;
  boolean cleaned_against_another_graph;
  StrBuf cleaned_against_graph_name;
} ErrorCleaning;

typedef struct
{
  StrBuf sample_name;
  uint64_t total_sequence;
  uint32_t mean_read_length;
  long double seq_err;
  ErrorCleaning cleaning;
} GraphInfo;

typedef struct
{
  HashTable ht;
  const uint32_t kmer_size, num_of_cols;
  uint32_t num_of_cols_used; // how many colours currently used

  GraphInfo *ginfo;

  // Optional fields:
  Edges *edges; // [hkey]

  // Colour specific arrays
  // cast to 2d array with:
  // Edges (*col_edges)[graph->num_of_cols]
  //   = (Edges (*)[graph->num_of_cols])graph->col_edges;
  // Covg (*col_covgs)[graph->num_of_cols]
  //   = (Covg (*)[graph->num_of_cols])graph->col_covgs;
  // then access with col_edges[hkey][col]
  Edges *col_edges; // [hkey*num_of_colours + col] or [hkey][col]
  Covg *col_covgs; // [hkey*num_of_colours + col] or [hkey][col]

  // [hkey/64][col] >> hkey%64
  // [num_of_colours*hkey/64+col] >> hkey%64
  uint64_t *node_in_cols;

  // path data
  uint64_t *kmer_paths;
  binary_paths_t pdata;

  // Loading reads
  uint64_t *readstrt;
} dBGraph;

#endif
