#ifndef CORRECT_ALN_INPUT_H_
#define CORRECT_ALN_INPUT_H_

#include "seqout.h"
#include "cortex_types.h"
#include "correct_alignment.h"
#include "async_read_io.h"

typedef struct
{
  AsyncIOInput files;
  uint8_t fq_cutoff, hp_cutoff;
  ReadMateDir matedir;
  CorrectAlnParam crt_params;
  // Next two only set if outputting sequences per file, as in ctx_correct.c
  char *out_base;
  SeqOutput *output;
} CorrectAlnInput;

#define CORRECT_ALN_INPUT_INIT {.fq_cutoff = 0, .hp_cutoff = 0,       \
                                .matedir = READPAIR_FR,               \
                                .crt_params = CORRECT_PARAMS_DEFAULT, \
                                .out_base = NULL, .output = NULL}

#include "objbuf_macro.h"
create_objbuf(correct_aln_input_buf, CorrectAlnInputBuffer, CorrectAlnInput);

void correct_aln_input_print(const CorrectAlnInput *c);

// Copy CorrectAlnInput to an array of AsyncIOInputs
void correct_aln_input_to_asycio(AsyncIOInput *asyncio_tasks,
                                 CorrectAlnInput *inputs,
                                 size_t num_inputs);

#endif /* CORRECT_ALN_INPUT_H_ */
