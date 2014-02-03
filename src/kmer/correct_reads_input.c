#include "global.h"
#include "correct_reads_input.h"

void correct_reads_input_init(const char *p1, const char *p2,
                              size_t col, size_t min_ins, size_t max_ins,
                              boolean one_way_bridge, size_t max_context,
                              float gap_variance, size_t gap_wiggle,
                              uint32_t fq_offset, uint32_t fq_cutoff,
                              uint32_t hp_cutoff, CorrectReadsInput *ptr)
{
  if(p1[0] == '-')
    die("Path appears to be an option: %s", p1);
  if(p2 != NULL && p2[0] == '-')
    die("Path appears to be an option: %s", p2);

  seq_file_t *f1, *f2 = NULL;

  if((f1 = seq_open(p1)) == NULL)
    die("Cannot read first --seq%s file: %s", p2 == NULL ? "" : "2", p1);
  if(p2 != NULL && (f2 = seq_open(p2)) == NULL)
    die("Cannot read second --seq2 file: %s", p2);

  CorrectReadsInput tsk = {.file1 = f1, .file2 = f2,
                           .ctxcol = 0, .ctpcol = col,
                           .ins_gap_min = min_ins, .ins_gap_max = max_ins,
                           .max_context = max_context,
                           .gap_variance = gap_variance,
                           .gap_wiggle = gap_wiggle,
                           .fq_offset = (uint8_t)fq_offset,
                           .fq_cutoff = (uint8_t)fq_cutoff,
                           .hp_cutoff = (uint8_t)hp_cutoff,
                           .read_pair_FR = true,
                           .one_way_gap_traverse = one_way_bridge};

  memcpy(ptr, &tsk, sizeof(CorrectReadsInput));
}

void correct_reads_input_print(const CorrectReadsInput *c)
{
  int has_p2 = c->file2 != NULL;
  const char *p1 = c->file1->path, *p2 = has_p2 ? c->file2->path : "";
  char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

  if(c->fq_cutoff > 0) sprintf(fqCutoff, "%u", c->fq_cutoff);
  if(c->fq_offset > 0) sprintf(fqOffset, "%u", c->fq_offset);
  if(c->hp_cutoff > 0) sprintf(hpCutoff, "%u", c->hp_cutoff);

  status("[task] input: %s%s%s", p1, (has_p2 ? ", " : ""), p2);
  status("  FASTQ offset: %s, threshold: %s; cut homopolymers: %s",
         fqOffset, fqCutoff, hpCutoff);

  // All one line
  timestamp();
  message("  %s-way gap traversal", c->one_way_gap_traverse ? "one" : "two");
  if(has_p2) {
    message("; read pair: %s; insert min,max: (%u,%u)",
            (c->read_pair_FR ? " FR" : "FF"), c->ins_gap_min, c->ins_gap_max);
  }
  message("\n");
}

// Sort by ctp colour, then by pointer address
static int correct_reads_input_cmp(const void *aa, const void *bb)
{
  const CorrectReadsInput *a = (const CorrectReadsInput*)aa;
  const CorrectReadsInput *b = (const CorrectReadsInput*)bb;
  if(a->ctpcol != b->ctpcol) return (int)a->ctpcol - (int)b->ctpcol;
  return (a > b ? 1 : (a < b ? -1 : 0));
}

void correct_reads_input_sort(CorrectReadsInput *inputs, size_t n)
{
  qsort(inputs, n, sizeof(CorrectReadsInput), correct_reads_input_cmp);
}

void correct_reads_input_to_asycio(AsyncIOReadTask *asyncio_tasks,
                                   CorrectReadsInput *inputs,
                                   size_t num_inputs)
{
  size_t i;
  for(i = 0; i < num_inputs; i++) {
    AsyncIOReadTask aio_task = {.file1 = inputs[i].file1,
                                .file2 = inputs[i].file2,
                                .fq_offset = inputs[i].fq_offset,
                                .ptr = &inputs[i]};

    memcpy(&asyncio_tasks[i], &aio_task, sizeof(AsyncIOReadTask));
  }
}
