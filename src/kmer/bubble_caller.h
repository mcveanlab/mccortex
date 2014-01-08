#ifndef BUBBLE_CALLER_H_
#define BUBBLE_CALLER_H_

#include "db_graph.h"
#include "cmd.h"

void bubble_caller_print_header(const dBGraph *db_graph, gzFile out,
                                const char* out_file, const CmdArgs *args);

// max_allele_len, max_flank_len in kmers
void invoke_bubble_caller(const dBGraph *db_graph, gzFile gzout,
                          size_t num_threads, char **tmp_paths,
                          size_t max_allele_len, size_t max_flank_len,
                          const size_t *ref_cols, size_t num_ref);

#endif /* BUBBLE_CALLER_H_ */
