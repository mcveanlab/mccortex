#include "global.h"

#include "seq_file.h"

#include "correct_reads_input.h"
#include "util.h"
#include "cmd.h"
#include "seq_reader.h"

// Should we print all paths in src/tools/generate_paths.c
bool gen_paths_print_contigs = false, gen_paths_print_paths = false;
bool gen_paths_print_reads = false;

void correct_reads_input_init(const char *p1, const char *p2, bool interleaved,
                              uint32_t fq_offset, uint32_t fq_cutoff,
                              uint32_t hp_cutoff, ReadMateDir matedir,
                              CorrectAlnParam params,
                              CorrectAlnReadsTask *ptr)
{
  // interleaved => p2 == NULL
  ctx_assert(!interleaved || p2 == NULL);

  if(p1[0] == '-')
    die("Path appears to be an option: %s", p1);
  if(p2 != NULL && p2[0] == '-')
    die("Path appears to be an option: %s", p2);

  seq_file_t *f1, *f2 = NULL;

  if((f1 = seq_open(p1)) == NULL)
    die("Cannot read first --seq%s file: %s", p2 == NULL ? "" : "2", p1);
  if(p2 != NULL && (f2 = seq_open(p2)) == NULL)
    die("Cannot read second --seq2 file: %s", p2);

  if(params.ins_gap_min > params.ins_gap_max) {
    die("--minIns %zu is greater than --maxIns %zu",
        (size_t)params.ins_gap_min, (size_t)params.ins_gap_max);
  }

  AsyncIOReadTask iotask = {.file1 = f1, .file2 = f2,
                            .fq_offset = (uint8_t)fq_offset,
                            .interleaved = interleaved,
                            .ptr = ptr};

  CorrectAlnReadsTask tsk = {.files = iotask,
                             .fq_cutoff = (uint8_t)fq_cutoff,
                             .hp_cutoff = (uint8_t)hp_cutoff,
                             .matedir = matedir,
                             .crt_params = params, .ptr = NULL};

  memcpy(ptr, &tsk, sizeof(CorrectAlnReadsTask));
}

void correct_reads_input_print(const CorrectAlnReadsTask *c)
{
  const AsyncIOReadTask *io = &c->files;
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

// Sort by ctp colour, then by pointer address
static int correct_reads_input_cmp(const void *aa, const void *bb)
{
  const CorrectAlnReadsTask *a = (const CorrectAlnReadsTask*)aa;
  const CorrectAlnReadsTask *b = (const CorrectAlnReadsTask*)bb;
  size_t col0 = a->crt_params.ctpcol, col1 = b->crt_params.ctpcol;
  if(col0 != col1) return (int)col0 - (int)col1;
  return (a > b ? 1 : (a < b ? -1 : 0));
}

void correct_reads_input_sort(CorrectAlnReadsTask *inputs, size_t n)
{
  qsort(inputs, n, sizeof(CorrectAlnReadsTask), correct_reads_input_cmp);
}

void correct_reads_input_to_asycio(AsyncIOReadTask *asyncio_tasks,
                                   CorrectAlnReadsTask *inputs,
                                   size_t num_inputs)
{
  size_t i;
  for(i = 0; i < num_inputs; i++) {
    memcpy(&asyncio_tasks[i], &inputs[i].files, sizeof(AsyncIOReadTask));
  }
}

#define DEFAULT_MIN_INS 0
#define DEFAULT_MAX_INS 500

// seq gap of N bases can be filled by MAX2(0, N±(N*GAP_VARIANCE+GAP_WIGGLE))
#define GAP_VARIANCE 0.1
#define GAP_WIGGLE 5

#define MAX_CONTEXT 200

// returns index of next argv
// use_pe means we are trying to bridge the insert gap, so take --(min|max)Ins
//
// without out_arg: --seq <in>       --seq2 <in1> <in2>       --seqi <in>
// with out_arg:    --seq <in> <out> --seq2 <in1> <in2> <out> --seqi <in> <out>
// Other options:
//  --fq_offset <int>
//  --fq_threshold <int>
//  --cut_hp <int>
//  --col <int>
//  --minIns <int>  if use_pe
//  --maxIns <int>  if use_pe
//  --seqgaps <out.csv> if dump_seq_gaps != NULL
//  --mpgaps <out.csv>  if dump_mp_gaps != NULL
//  --oneway
//  --twoway
//  --FF --FR --RF --RR
//  --no-end-check option // --end-check
int correct_reads_parse(int argc, char **argv,
                        bool use_pe, bool out_arg,
                        CorrectAlnReadsTask *inputs,
                        size_t *num_inputs_ptr,
                        char **dump_seq_gaps, char **dump_mp_gaps)
{
  int argi;
  size_t num_inputs = 0, col, min_ins, max_ins;
  uint32_t fq_offset = 0, fq_cutoff = 0, hp_cutoff = 0;
  bool col_set = false, col_used = false;
  ReadMateDir matedir = READPAIR_FR;

  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .ins_gap_min = DEFAULT_MIN_INS,
                            .ins_gap_max = DEFAULT_MAX_INS,
                            .one_way_gap_traverse = true,
                            .use_end_check = true,
                            .max_context = MAX_CONTEXT,
                            .gap_variance = GAP_VARIANCE,
                            .gap_wiggle = GAP_WIGGLE};

  // These just used for printing warnings / checking arg combinations
  bool seen_pe_reads = false, seen_dump_mp_gaps = false;

  for(argi = 0; argi < argc && argv[argi][0] == '-' &&argv[argi][1]; argi++)
  {
    if(!strcmp(argv[argi],"--fq_offset") || !strcmp(argv[argi],"-q")) {
      if(argi + 1 >= argc)
        cmd_print_usage("--fq_offset <offset> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_offset) || fq_offset > 128)
        die("Invalid --fq_offset argument: %s", argv[argi+1]);
      argi++;
    }
    else if(!strcmp(argv[argi],"--fq_threshold") || !strcmp(argv[argi],"-Q")) {
      if(argi + 1 >= argc)
        cmd_print_usage("--fq_threshold <qual> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_offset) || fq_offset > 128)
        die("Invalid --fq_threshold argument: %s", argv[argi+1]);
      argi++;
    }
    else if(!strcmp(argv[argi],"--cut_hp") || !strcmp(argv[argi],"-h")) {
      if(argi + 1 >= argc)
        cmd_print_usage("--cut_hp <len> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &hp_cutoff))
        die("Invalid --cut_hp argument: %s", argv[argi+1]);
      if(hp_cutoff > UINT8_MAX)
        die("--cut_hp <hp> cannot be greater than %i", UINT8_MAX);
      argi++;
    }
    else if(!strcmp(argv[argi],"--col") || !strcmp(argv[argi],"--colour") ||
            !strcmp(argv[argi],"--color") || !strcmp(argv[argi],"-c"))
    {
      if(argi+2 >= argc) cmd_print_usage("--col <colour> requires an argument");
      if(col_set && !col_used)
        cmd_print_usage("--seq or --seq2 must follow --col");
      if(!parse_entire_size(argv[argi+1], &col))
        cmd_print_usage("--col <colour> requires integers >= 0");
      params.ctpcol = params.ctxcol = col;
      col_set = true;
      col_used = false;
      argi++;
    }
    else if(use_pe && (!strcmp(argv[argi],"--minIns") || !strcmp(argv[argi], "-i")))
    {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &min_ins))
        cmd_print_usage("--minIns <bp> requires a positive integer arg");
      params.ins_gap_min = min_ins;
      argi++;
    }
    else if(use_pe && (!strcmp(argv[argi],"--maxIns") || !strcmp(argv[argi], "-j")))
    {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &max_ins))
        cmd_print_usage("--maxIns <bp> requires a positive integer arg");
      if(max_ins < 20)
        warn("--maxGap < 20 seems very low!");
      params.ins_gap_max = max_ins;
      argi++;
    }
    else if(!strcmp(argv[argi],"--end-check") || !strcmp(argv[argi],"-e"))
    {
      params.use_end_check = true;
      argi++;
    }
    else if(!strcmp(argv[argi],"--no-end-check") || !strcmp(argv[argi],"-E"))
    {
      params.use_end_check = false;
      argi++;
    }
    else if(dump_seq_gaps && (!strcmp(argv[argi],"--seqgaps") || !strcmp(argv[argi],"-g")))
    {
      if(argi+1 >= argc) cmd_print_usage("--seqgaps <out.csv> requires arg");
      *dump_seq_gaps = argv[argi+1];
      argi++;
    }
    else if(dump_mp_gaps && (!strcmp(argv[argi],"--mpgaps") || !strcmp(argv[argi],"-G")))
    {
      if(argi+1 >= argc) cmd_print_usage("--mpgaps <out.csv> requires arg");
      *dump_mp_gaps = argv[argi+1];
      argi++;
      seen_dump_mp_gaps = true;
    }
    else if(!strcmp(argv[argi],"--FF")) matedir = READPAIR_FF;
    else if(!strcmp(argv[argi],"--FR")) matedir = READPAIR_FR;
    else if(!strcmp(argv[argi],"--RF")) matedir = READPAIR_RF;
    else if(!strcmp(argv[argi],"--RR")) matedir = READPAIR_RR;
    else if(!strcmp(argv[argi],"--oneway") || !strcmp(argv[argi],"-w"))
      params.one_way_gap_traverse = true;
    else if(!strcmp(argv[argi],"--twoway") || !strcmp(argv[argi],"-W"))
      params.one_way_gap_traverse = false;
    else if(!strcmp(argv[argi],"--seq") || !strcmp(argv[argi],"-1") ||
            !strcmp(argv[argi],"--seq") || !strcmp(argv[argi],"-i"))
    {
      if(out_arg) {
        if(argi+2 >= argc) cmd_print_usage("--seq <in> <out> missing args");
        if(!col_set) die("--seq <in.fa> <out> before --col <colour>");
      } else {
        if(argi+1 >= argc) cmd_print_usage("--seq <in> missing args");
        if(!col_set) die("--seq <in.fa> before --col <colour>");
      }

      bool interleaved = (!strcmp(argv[argi],"--seqi"));

      correct_reads_input_init(argv[argi+1], NULL, interleaved,
                               fq_offset, fq_cutoff, hp_cutoff, matedir, params,
                               &inputs[num_inputs]);

      if(out_arg) {
        // Use void pointer to store output path
        inputs[num_inputs].ptr = argv[argi+2];
      }

      num_inputs++;
      col_used = true;
      argi += 1+out_arg;
    }
    else if(!strcmp(argv[argi],"--seq2") || !strcmp(argv[argi],"-2"))
    {
      if(out_arg) {
        if(argi+3 >= argc) cmd_print_usage("--seq2 <in1> <in2> <out> missing args");
        if(!col_set) die("--seq2 <in1> <in2> <out> before --col <colour>");
      } else {
        if(argi+2 >= argc) cmd_print_usage("--seq2 <in1> <in2> missing args");
        if(!col_set) die("--seq2 <in1> <in2> before --col <colour>");
      }

      correct_reads_input_init(argv[argi+1], argv[argi+2], false,
                               fq_offset, fq_cutoff, hp_cutoff, matedir, params,
                               &inputs[num_inputs]);

      if(out_arg) {
        // Use void pointer to store output path
        inputs[num_inputs].ptr = argv[argi+3];
      }

      num_inputs++;
      col_used = true;
      argi += 2+out_arg;

      seen_pe_reads = true;
    }
    else cmd_print_usage("Unknown option: %s", argv[argi]);
  }

  if(num_inputs == 0)
    cmd_print_usage("need at least one: --col <c> [--seq|--seq2]");
  if(!col_used)
    cmd_print_usage("--seq must follow last --col <col>");

  correct_reads_input_sort(inputs, num_inputs);

  if(!seen_pe_reads && seen_dump_mp_gaps)
    warn("--mpgaps without paired end reads from --seq2");

  *num_inputs_ptr = num_inputs;
  return argi;
}
