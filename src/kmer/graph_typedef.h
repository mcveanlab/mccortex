#ifndef GRAPH_TYPEDEF_H_
#define GRAPH_TYPEDEF_H_

#include "string_buffer.h"

//
// This file holds type definitions used in db_graph
//

typedef size_t Colour;
typedef uint8_t Edges;
typedef uint32_t Covg;

#define COVG_MAX UINT_MAX

#define FORWARD 0
#define REVERSE 1
typedef uint8_t Orientation;

//
// Paths
//

typedef uint64_t PathIndex;
typedef uint16_t PathLen;

typedef struct {
  uint8_t *const store, *const end;
  const size_t size;
  const size_t num_of_cols, colset_bytes;
  uint8_t *next;
  size_t num_of_paths, num_kmers_with_paths;
} PathStore;

//
// Graph Info including cleaning
//

// Thesholds are zero if not used (e.g. remv_low_cov_sups == false)
// is_graph_intersection is for cleaning a low covg sample against
// cleaned pool of population
typedef struct
{
  boolean tip_clipping, remv_low_cov_sups, remv_low_cov_nodes;
  Covg remv_low_cov_sups_thresh, remv_low_cov_nodes_thresh;
  boolean is_graph_intersection; // formerly cleaned_against_another_graph
  StrBuf intersection_name; // formerly cleaned_against_graph_nme
} ErrorCleaning;

typedef struct
{
  uint32_t mean_read_length;
  uint64_t total_sequence;
  StrBuf sample_name;
  long double seq_err;
  ErrorCleaning cleaning;
} GraphInfo;

#endif
