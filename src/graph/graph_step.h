#ifndef GRAPH_STEP_H_
#define GRAPH_STEP_H_

/*
  This file contains the struct and constants used to record the behaviour of
  the GraphWalker at each "step".
*/

// GraphStep.status values:
#define GRPHWLK_FORWARD       0 /* Success: only one choice */
#define GRPHWLK_COLFWD        1 /* Success: only one choice in colour */
#define GRPHWLK_NOCOVG        2 /* Fail: no choices */
#define GRPHWLK_NOCOLCOVG     3 /* Fail: fork in pop but no choices in colour */
#define GRPHWLK_NOPATHS       4 /* Fail: fork in colour, no paths */
#define GRPHWLK_SPLIT_PATHS   5 /* Fail: fork in colour, paths split */
#define GRPHWLK_MISSING_PATHS 6 /* Fail: fork in colour, missing info */
#define GRPHWLK_USEPATH       7 /* Success: fork in colour, paths resolved */
#define GRPHWLK_NUM_STATES    8

#define GRPHWLK_FORWARD_STR       "GoForward"
#define GRPHWLK_COLFWD_STR        "GoColForward"
#define GRPHWLK_NOCOVG_STR        "FailNoCovg"
#define GRPHWLK_NOCOLCOVG_STR     "FailNoColCovg"
#define GRPHWLK_NOPATHS_STR       "FailNoPaths"
#define GRPHWLK_SPLIT_PATHS_STR   "FailSplitPaths"
#define GRPHWLK_MISSING_PATHS_STR "FailMissingPaths"
#define GRPHWLK_USEPATH_STR       "GoUsePath"

// Result from graph_walker_choose
typedef struct
{
  // idx is -1 if failed, otherwise index of node taken [0..3]
  int8_t idx;
  uint8_t status;
  bool node_has_col;
  // path_step is number of nodes in the repeat that have been resolved
  // immediately after a junction. In order to resolve the repeat,
  // read required to resolve = path_gap+kmer_size-1+2 bp
  //   path_gap+kmer_size-1 is the number of bases in the repeat
  //   +2 if for a base either side to link a node either side of the tangle
  size_t path_gap;
} GraphStep;

extern const char *graph_step_str[];

// Was the last step resolving a split in this colour?
#define graph_step_status_is_fork(stat) ((stat) > GRPHWLK_NOCOLCOVG)

// Are we still walking?
#define grap_step_status_is_good(stat) ((stat) <= GRPHWLK_COLFWD || \
                                        (stat) == GRPHWLK_USEPATH)

char* graph_step_status2str(uint8_t status, char *str, size_t len);
void graph_step_print_state_hist(const size_t arr[GRPHWLK_NUM_STATES]);

#endif /* GRAPH_STEP_H_ */
