#include "global.h"
#include "correct_aln_input.h"
#include "seq_reader.h"
#include "file_util.h"

cJSON* correct_aln_input_json_hdr(const CorrectAlnInput *input)
{
  cJSON *hdr = cJSON_CreateObject();

  const char *files[2] = {input->files.file1->path, NULL};
  int nfiles = 1;
  if(input->files.file2) { files[1] = input->files.file2->path; nfiles++; }
  cJSON_AddItemToObject(hdr, "files", cJSON_CreateStringArray(files, nfiles));

  cJSON_AddItemToObject(hdr, "interleaved", cJSON_CreateBool(input->files.interleaved));
  cJSON_AddItemToObject(hdr, "fq_offset", cJSON_CreateInt(input->files.fq_offset));
  cJSON_AddItemToObject(hdr, "fq_cutoff", cJSON_CreateInt(input->fq_cutoff));
  cJSON_AddItemToObject(hdr, "hp_cutoff", cJSON_CreateInt(input->hp_cutoff));
  cJSON_AddItemToObject(hdr, "matepair", cJSON_CreateString(MP_DIR_STRS[input->matedir]));

  const CorrectAlnParam *p = &input->crt_params;
  cJSON_AddItemToObject(hdr, "frag_len_min_bp",  cJSON_CreateInt(p->frag_len_min));
  cJSON_AddItemToObject(hdr, "frag_len_max_bp",  cJSON_CreateInt(p->frag_len_max));
  cJSON_AddItemToObject(hdr, "one_way_gap_fill", cJSON_CreateBool(p->one_way_gap_traverse));
  cJSON_AddItemToObject(hdr, "use_end_check",    cJSON_CreateBool(p->use_end_check));
  cJSON_AddItemToObject(hdr, "max_context",      cJSON_CreateInt(p->max_context));
  cJSON_AddItemToObject(hdr, "gap_variance",     cJSON_CreateNumber(p->gap_variance));
  cJSON_AddItemToObject(hdr, "gap_wiggle",       cJSON_CreateNumber(p->gap_wiggle));

  return hdr;
}

void correct_aln_input_print(const CorrectAlnInput *c)
{
  const AsyncIOInput *io = &c->files;
  int has_p2 = io->file2 != NULL;
  const char *p1 = io->file1->path, *p2 = has_p2 ? io->file2->path : "";
  char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

  if(io->fq_offset > 0) sprintf(fqOffset, "%u", io->fq_offset);
  if(c->fq_cutoff > 0) sprintf(fqCutoff, "%u", c->fq_cutoff);
  if(c->hp_cutoff > 0) sprintf(hpCutoff, "%u", c->hp_cutoff);

  status("[task] input: %s%s%s", futil_inpath_str(p1),
         (has_p2 ? ", " : ""), futil_inpath_str(p2));
  status("  FASTQ offset: %s, threshold: %s; cut homopolymers: %s",
         fqOffset, fqCutoff, hpCutoff);

  // All one line
  timestamp();
  message("  %s-way gap traversal", c->crt_params.one_way_gap_traverse ? "one" : "two");
  if(has_p2) {
    message("; read pair: %s; insert min,max: (%u,%u)", MP_DIR_STRS[c->matedir],
            c->crt_params.frag_len_min, c->crt_params.frag_len_max);
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
