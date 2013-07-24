#ifndef GRAPH_TYPEDEF_H_
#define GRAPH_TYPEDEF_H_

#include "string_buffer.h"
#include "hash_table.h"

typedef uint8_t Colour;
typedef uint8_t Edges;
typedef uint32_t Covg;

#define COVG_MAX UINT_MAX

#define FORWARD 0
#define REVERSE 1
typedef uint8_t Orientation;

typedef BinaryKmerPtr Key;

// Paths

#define PATH_NULL UINT64_MAX
#define PATH_LEN_BITS 15
#define MAX_PATHLEN ((1U<<PATH_LEN_BITS)-1)

typedef uint64_t PathIndex;
typedef uint16_t PathLen;

typedef struct {
  uint8_t *const store, *const end;
  const size_t size;
  const uint32_t num_of_cols, col_bitset_bytes;
  uint8_t *next;
  size_t num_of_paths, num_kmers_with_paths;
} PathStore;

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
  uint32_t mean_read_length;
  uint64_t total_sequence;
  StrBuf sample_name;
  long double seq_err;
  ErrorCleaning cleaning;
} GraphInfo;

typedef struct
{
  HashTable ht;
  const uint32_t kmer_size, num_of_cols;
  uint32_t num_of_cols_used; // how many colours currently used

  // Array of GraphInfo objects, one per colour
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
  volatile PathIndex *kmer_paths;
  PathStore pdata;

  // Loading reads
  uint64_t *readstrt;
} dBGraph;

#endif
