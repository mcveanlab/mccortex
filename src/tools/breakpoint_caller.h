#ifndef BREAKPOINT_CALLER_H_
#define BREAKPOINT_CALLER_H_

#include "seq_file.h"
#include "db_graph.h"
#include "cmd.h"

#include "cJSON/cJSON.h"

#define BREAKPOINT_FORMAT_VERSION 3

// Require 5 kmers on the reference before and after breakpoint
#define DEFAULT_MIN_REF_NKMERS 5

/**
 * Adds input bkmers to the graph
 * @param hdrs JSON headers of input files
 * @param ref_col colour to add reference kmers to
 **/
void breakpoints_call(size_t num_of_threads, size_t ref_col,
                      gzFile gzout, const char *out_path,
                      const read_t *reads, size_t num_reads,
                      char **seq_paths, size_t num_seq_paths,
                      size_t min_ref_flank,
                      cJSON **hdrs, size_t nhdrs,
                      dBGraph *db_graph);

#endif /* BREAKPOINT_CALLER_H_ */
