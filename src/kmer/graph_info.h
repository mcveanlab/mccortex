#ifndef DB_INFO_H_
#define DB_INFO_H_

#include <inttypes.h>
#include "string_buffer.h"
#include "cortex_types.h"

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

void graph_info_init(GraphInfo *ginfo);
void graph_info_alloc(GraphInfo *ginfo);
void graph_info_dealloc(GraphInfo *ginfo);

void graph_info_make_intersect(const GraphInfo *ginfo, StrBuf *intersect_name);
void graph_info_append_intersect(ErrorCleaning *cleaning, const char *intersect_name);

// Update mean read length in a colour, eg when you merge a new binary
// what it used to be, ie the previous_seq)
void graph_info_update_contigs(GraphInfo *ginfo,
                               uint64_t added_seq, uint64_t num_contigs);

void graph_info_cpy(GraphInfo *dst, const GraphInfo *src);
void graph_info_merge(GraphInfo *dst, const GraphInfo *src);

#endif /* GRAPH_INFO_H_ */
