#ifndef CORRECT_READS_H_
#define CORRECT_READS_H_

#include "correct_aln_input.h"
#include "db_graph.h"

// Correct reads against the graph, and print out
void correct_reads(CorrectAlnInput *inputs, size_t num_inputs,
                   const char *dump_seqgap_hist_path,
                   const char *dump_fraglen_hist_path,
                   size_t num_threads, const dBGraph *db_graph);

#endif /* CORRECT_READS_H_ */
