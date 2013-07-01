#ifndef DB_INFO_H_
#define DB_INFO_H_

#include <inttypes.h>
#include "graph_typedef.h"

void error_cleaning_init(ErrorCleaning *ec);
void error_cleaning_alloc(ErrorCleaning *ec);
void error_cleaning_dealloc(ErrorCleaning *ec);

void error_cleaning_overwrite(ErrorCleaning *tgt, const ErrorCleaning *src);

void graph_info_init(GraphInfo *ginfo);
void graph_info_alloc(GraphInfo *ginfo);
void graph_info_dealloc(GraphInfo *ginfo);

// Update mean read length in a colour, eg when you merge a new binary
// what it used to be, ie the previous_seq)
void graph_info_update_seq_stats(GraphInfo *ginfo,
                                 uint32_t added_mean, uint64_t added_seq);

#endif /* GRAPH_INFO_H_ */
