#include "global.h"
#include "graph_step.h"
#include "util.h"

/*
  This file contains the struct and constants used to record the behaviour of
  the GraphWalker at each "step".
*/

const char *graph_step_str[] = {GRPHWLK_FORWARD_STR,       GRPHWLK_COLFWD_STR,
                                GRPHWLK_NOCOVG_STR,        GRPHWLK_NOCOLCOVG_STR,
                                GRPHWLK_NOPATHS_STR,       GRPHWLK_SPLIT_PATHS_STR,
                                GRPHWLK_MISSING_PATHS_STR, GRPHWLK_USEPATH_STR};

char* graph_step_status2str(uint8_t status, char *str, size_t len)
{
  ctx_assert(len >= 20); (void)len;
  ctx_assert(status < GRPHWLK_NUM_STATES);
  strcpy(str, graph_step_str[status]);
  return str;
}

void graph_step_print_state_hist(const size_t hist[GRPHWLK_NUM_STATES])
{
  util_print_nums(graph_step_str, hist, GRPHWLK_NUM_STATES, 30);
}
