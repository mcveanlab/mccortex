#ifndef BREAKPOINT_CALLER_H_
#define BREAKPOINT_CALLER_H_

#include "seq_file.h"
#include "db_graph.h"
#include "cmd.h"

// Adds read bkmers to the graph
void breakpoints_call(gzFile gzout, dBGraph *db_graph,
                      const read_t *reads, size_t num_reads,
                      size_t num_of_threads, const CmdArgs *args);

#endif /* BREAKPOINT_CALLER_H_ */
