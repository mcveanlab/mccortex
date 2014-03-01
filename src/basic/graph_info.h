#ifndef DB_INFO_H_
#define DB_INFO_H_

#include <inttypes.h>
#include "string_buffer.h"
#include "cortex_types.h"
#include "loading_stats.h"

// Thesholds are zero if not used (e.g. cleaned_snodes == false)
// is_graph_intersection is for cleaning a low covg sample against
// cleaned pool of population
typedef struct
{
  bool cleaned_tips, cleaned_snodes, cleaned_kmers;
  Covg clean_snodes_thresh, clean_kmers_thresh;
  bool is_graph_intersection;
  StrBuf intersection_name;
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
void graph_info_append_intersect(ErrorCleaning *cleaning,
                                 const char *intersect_name);

void graph_info_cpy(GraphInfo *dst, const GraphInfo *src);
void graph_info_merge(GraphInfo *dst, const GraphInfo *src);

void graph_info_update_stats(GraphInfo *ginfo, const LoadingStats *stats);

#endif /* GRAPH_INFO_H_ */
