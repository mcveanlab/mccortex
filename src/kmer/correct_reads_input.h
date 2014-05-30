#ifndef CORRECT_READ_INPUT_H_
#define CORRECT_READ_INPUT_H_

#include "seq_file.h"
#include "cortex_types.h"
#include "correct_alignment.h"
#include "async_read_io.h"

typedef struct
{
  AsyncIOReadTask files;
  uint8_t fq_cutoff, hp_cutoff;
  ReadMateDir matedir;
  CorrectAlnParam crt_params;
  void *ptr; // general porpoise pointer
} CorrectAlnTask;


#include "objbuf_macro.h"
create_objbuf(correct_aln_task_buf, CorrectAlnTaskBuffer, CorrectAlnTask);

// These are used in generate_paths.c
extern bool gen_paths_print_contigs, gen_paths_print_paths, gen_paths_print_reads;

int correct_reads_parse(int argc, char **argv,
                        bool use_pe, bool out_arg,
                        CorrectAlnTask *inputs, size_t *num_inputs_ptr,
                        char **dump_seq_gaps, char **dump_mp_gaps);

void correct_reads_input_print(const CorrectAlnTask *c);

void correct_reads_input_sort(CorrectAlnTask *inputs, size_t n);

// Copy CorrectAlnTask to an array of AsyncIOReadTasks
void correct_reads_input_to_asycio(AsyncIOReadTask *asyncio_tasks,
                                   CorrectAlnTask *inputs,
                                   size_t num_inputs);

#endif /* CORRECT_READ_INPUT_H_ */
