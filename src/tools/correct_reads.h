#ifndef CORRECT_READS_H_
#define CORRECT_READS_H_

#include "correct_aln_input.h"
#include "db_graph.h"

// Correct reads against the graph, and print out
// `input` and `outputs` should both be of length `num_inputs`
void correct_reads(size_t num_threads, size_t max_io_threads,
                   CorrectAlnInput *inputs, size_t num_inputs,
                   const dBGraph *db_graph);

#endif /* CORRECT_READS_H_ */
