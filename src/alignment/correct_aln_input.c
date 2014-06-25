#include "global.h"
#include "correct_aln_input.h"
#include "seq_reader.h"
#include "file_util.h"

void correct_aln_input_print(const CorrectAlnInput *c)
{
  const AsyncIOInput *io = &c->files;
  int has_p2 = io->file2 != NULL;
  const char *p1 = io->file1->path, *p2 = has_p2 ? io->file2->path : "";
  char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

  if(io->fq_offset > 0) sprintf(fqOffset, "%u", io->fq_offset);
  if(c->fq_cutoff > 0) sprintf(fqCutoff, "%u", c->fq_cutoff);
  if(c->hp_cutoff > 0) sprintf(hpCutoff, "%u", c->hp_cutoff);

  status("[task] input: %s%s%s", p1, (has_p2 ? ", " : ""), p2);
  status("  FASTQ offset: %s, threshold: %s; cut homopolymers: %s",
         fqOffset, fqCutoff, hpCutoff);

  // All one line
  timestamp();
  message("  %s-way gap traversal", c->crt_params.one_way_gap_traverse ? "one" : "two");
  if(has_p2) {
    message("; read pair: %s; insert min,max: (%u,%u)", MP_DIR_STRS[c->matedir],
            c->crt_params.ins_gap_min, c->crt_params.ins_gap_max);
  }
  message(" [%sedge check]", c->crt_params.use_end_check ? "" : "no ");
  message("\n");
}

void correct_aln_input_to_asycio(AsyncIOInput *asyncio_tasks,
                                 CorrectAlnInput *inputs,
                                 size_t num_inputs)
{
  size_t i;
  for(i = 0; i < num_inputs; i++) {
    memcpy(&asyncio_tasks[i], &inputs[i].files, sizeof(AsyncIOInput));
  }
}
