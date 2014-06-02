#ifndef CORRECT_READS_H_
#define CORRECT_READS_H_

#include "correct_reads_input.h"
#include "db_graph.h"

typedef struct {
  StrBuf outse, outpe1, outpe2;
  gzFile gzse, gzpe1, gzpe2;
  pthread_mutex_t lockse, lockpe;
} CorrectedOutput;

// Returns true on success, false on error
bool corrected_output_open(CorrectedOutput *out, const CorrectAlnTask *in);
void corrected_output_close(CorrectedOutput *out);
void corrected_output_delete(CorrectedOutput *out);

// Correct reads against the graph, and print out
// `input` and `outputs` should both be of length `num_inputs`
void correct_reads(size_t num_threads, size_t max_io_threads,
                   CorrectAlnTask *inputs, CorrectedOutput *outputs,
                   size_t num_inputs, const dBGraph *db_graph);

#endif /* CORRECT_READS_H_ */
