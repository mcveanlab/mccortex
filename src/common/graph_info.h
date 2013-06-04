#ifndef DB_INFO_H_
#define DB_INFO_H_

#include <inttypes.h>
#include "string_buffer.h"

typedef struct
{
  boolean tip_clipping, remv_low_cov_sups, remv_low_cov_nodes;

  // Thesholds are zero if not used (e.g. remv_low_cov_sups == false)
  Covg remv_low_cov_sups_thresh, remv_low_cov_nodes_thresh;

  // Cleaning a low covg sample against cleaned pool of population
  boolean cleaned_against_another_graph;
  StrBuf* cleaned_against_graph_name;
} ErrorCleaning;

ErrorCleaning* error_cleaning_alloc(ErrorCleaning *ec);
void error_cleaning_dealloc(ErrorCleaning *ec);
void error_cleaning_init(ErrorCleaning *ec);

typedef struct
{
  StrBuf *sample_names[NUM_OF_COLOURS];
  uint64_t total_sequence[NUM_OF_COLOURS];
  uint32_t mean_read_length[NUM_OF_COLOURS];
  long double seq_err[NUM_OF_COLOURS];
  ErrorCleaning cleaning[NUM_OF_COLOURS];
  uint32_t num_of_colours_loaded, num_of_shades_loaded;
} GraphInfo;

GraphInfo* graph_info_alloc(GraphInfo *ginfo);
void graph_info_dealloc(GraphInfo *ginfo);
void graph_info_init(GraphInfo *ginfo);

void graph_info_reset_one_colour(GraphInfo *ginfo, uint32_t colour);

// Set total amount of sequence in a colour
void graph_info_set_seq(GraphInfo *ginfo, int colour, uint64_t num_bp);
void graph_info_set_mean_readlen(GraphInfo *ginfo, int colour, int len);
void graph_info_set_tip_clipping(GraphInfo *ginfo, int colour);
void graph_info_set_remv_low_cov_sups(GraphInfo *ginfo, int colour, Covg thresh);
void graph_info_set_remv_low_cov_nodes(GraphInfo *ginfo, int colour, Covg thresh);
void graph_info_set_seq_err(GraphInfo *ginfo, int col, long double err);

// Update mean read length in a colour, eg when you merge a new binary
// what it used to be, ie the previous_seq)
void graph_info_update_mean_readlen_and_total_seq(GraphInfo *ginfo,
                                                  uint32_t colour,
                                                  uint32_t added_mean,
                                                  uint64_t added_seq);

int get_mean_readlen_across_colours(GraphInfo *ginfo);
void read_estimated_seq_errors_from_file(GraphInfo *ginfo, FILE *fp);
void print_seq_err_rates_to_screen(GraphInfo *ginfo);

#endif /* GRAPH_INFO_H_ */
