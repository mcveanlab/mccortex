#ifndef CORRECT_READ_INPUT_H_
#define CORRECT_READ_INPUT_H_

#include "seq_file.h"
#include "cortex_types.h"
#include "async_read_io.h"

typedef struct
{
  seq_file_t *const file1, *const file2;
  const Colour ctpcol, ctxcol;
  const uint32_t ins_gap_min, ins_gap_max;
  const boolean one_way_gap_traverse; // set to false for more error prone algo
  const uint32_t max_context, gap_wiggle;
  const float gap_variance; // permitted gap size = X*gap_variance + gap_wiggle
  const uint8_t fq_offset, fq_cutoff, hp_cutoff;
  const boolean read_pair_FR; // set to false if reads are FF
  void *ptr;
} CorrectReadsInput;

void correct_reads_input_init(const char *p1, const char *p2,
                              size_t col, size_t min_ins, size_t max_ins,
                              boolean one_way_bridge, size_t max_context,
                              float gap_variance, size_t gap_wiggle,
                              uint32_t fq_offset, uint32_t fq_cutoff,
                              uint32_t hp_cutoff, CorrectReadsInput *ptr);

void correct_reads_input_print(const CorrectReadsInput *c);

void correct_reads_input_sort(CorrectReadsInput *inputs, size_t n);

// Copy CorrectReadsInput to an array of AsyncIOReadTasks
void correct_reads_input_to_asycio(AsyncIOReadTask *asyncio_tasks,
                                   CorrectReadsInput *inputs,
                                   size_t num_inputs);

#endif /* CORRECT_READ_INPUT_H_ */
