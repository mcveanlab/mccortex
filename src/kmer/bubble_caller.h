#ifndef BUBBLE_CALLER_H_
#define BUBBLE_CALLER_H_

#include "db_graph.h"
#include "cmd.h"

void invoke_bubble_caller(const dBGraph *db_graph, const char* out_file,
                          int num_threads, char **tmp_paths,
                          size_t max_allele_len, size_t max_flank_len,
                          const size_t *ref_cols, size_t num_ref,
                          const CmdArgs *args);

#endif /* BUBBLE_CALLER_H_ */
