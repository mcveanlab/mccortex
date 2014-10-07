#include "global.h"
#include "cmd.h"
#include "read_thread_cmd.h"
#include "gpath_checks.h"
#include "file_util.h"

//
// ctx_thread.c and ctx_correct.c use many of the same command line arguments
// Therefore this file provides the common functionality to both
//

void read_thread_args_alloc(struct ReadThreadCmdArgs *args)
{
  memset(args, 0, sizeof(*args));
  correct_aln_input_buf_alloc(&args->inputs, 16);
  gpfile_buf_alloc(&args->gpfiles, 16);
}

void read_thread_args_dealloc(struct ReadThreadCmdArgs *args)
{
  size_t i;
  for(i = 0; i < args->inputs.len; i++) asyncio_task_close(&args->inputs.data[i].files);
  for(i = 0; i < args->gpfiles.len; i++) gpath_reader_close(&args->gpfiles.data[i]);

  correct_aln_input_buf_dealloc(&args->inputs);
  gpfile_buf_dealloc(&args->gpfiles);
}

void read_thread_args_parse(struct ReadThreadCmdArgs *args,
                            int argc, char **argv,
                            const struct option *longopts, bool correct_cmd)
{
  size_t i;
  CorrectAlnInput task = CORRECT_ALN_INPUT_INIT;
  uint8_t fq_offset = 0;
  size_t dump_seq_hist_n = 0, dump_frag_hist_n = 0; // how many times are -g -G specified
  GPathReader tmp_gpfile;

  CorrectAlnInputBuffer *inputs = &args->inputs;
  args->memargs = (struct MemArgs)MEM_ARGS_INIT;
  args->fmt = SEQ_FMT_FASTQ;

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int used = 1, c;
  char *tmp_path;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!args->out_ctp_path,cmd); args->out_ctp_path = optarg; break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_add(&args->gpfiles, tmp_gpfile);
        break;
      case 't':
        cmd_check(!args->nthreads, cmd);
        args->nthreads = cmd_uint32_nonzero(cmd, optarg);
        break;
      case 'm': cmd_mem_args_set_memory(&args->memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&args->memargs, optarg); break;
      case 'c': args->colour = cmd_uint32(cmd, optarg); break;
      case 'F':
        cmd_check(args->fmt == SEQ_FMT_FASTQ, cmd);
        args->fmt = cmd_parse_format(cmd, optarg);
        break;
      case '1':
      case '2':
      case 'i':
        used = 1;
        correct_aln_input_buf_add(inputs, task);
        asyncio_task_parse(&inputs->data[inputs->len-1].files, c, optarg,
                           fq_offset, correct_cmd ? &tmp_path : NULL);
        if(correct_cmd) inputs->data[inputs->len-1].out_base = tmp_path;
        break;
      case 'M':
             if(!strcmp(optarg,"FF")) task.matedir = READPAIR_FF;
        else if(!strcmp(optarg,"FR")) task.matedir = READPAIR_FR;
        else if(!strcmp(optarg,"RF")) task.matedir = READPAIR_RF;
        else if(!strcmp(optarg,"RR")) task.matedir = READPAIR_RR;
        else die("-M,--matepair <orient> must be one of: FF,FR,RF,RR");
        used = 0; break;
      case 'O': fq_offset = cmd_uint8(cmd, optarg); used = 0; break;
      case 'Q': task.fq_cutoff = cmd_uint8(cmd, optarg); used = 0; break;
      case 'H': task.hp_cutoff = cmd_uint8(cmd, optarg); used = 0; break;
      case 'l': task.crt_params.frag_len_min = cmd_uint32(cmd, optarg); used = 0; break;
      case 'L': task.crt_params.frag_len_max = cmd_uint32(cmd, optarg); used = 0; break;
      case 'w': task.crt_params.one_way_gap_traverse = true; used = 0; break;
      case 'W': task.crt_params.one_way_gap_traverse = false; used = 0; break;
      case 'd': task.crt_params.gap_wiggle = cmd_udouble(cmd, optarg); used = 0; break;
      case 'D': task.crt_params.gap_variance = cmd_udouble(cmd, optarg); used = 0; break;
      case 'X': task.crt_params.max_context = cmd_uint32(cmd, optarg); used = 0; break;
      case 'e': task.crt_params.use_end_check = true; used = 0; break;
      case 'E': task.crt_params.use_end_check = false; used = 0; break;
      case 'g': args->dump_seq_sizes = optarg; dump_seq_hist_n++; break;
      case 'G': args->dump_frag_sizes = optarg; dump_frag_hist_n++; break;
      case 'u': args->use_new_paths = true; break;
      case 'x': gen_paths_print_contigs = true; break;
      case 'y': gen_paths_print_paths = true; break;
      case 'z': gen_paths_print_reads = true; break;
      case 'Z':
        cmd_check(!args->fq_zero, cmd);
        if(strlen(optarg) != 1)
          cmd_print_usage("--fq-zero <c> requires a single char");
        args->fq_zero = optarg[0];
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" thread/correct -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(args->nthreads == 0) args->nthreads = DEFAULT_NTHREADS;

  // Check that optind+1 == argc
  if(optind+1 > argc)
    cmd_print_usage("Expected exactly one graph file");
  else if(optind+1 < argc)
    cmd_print_usage("Expected only one graph file. What is this: '%s'", argv[optind]);

  char *graph_path = argv[optind];
  status("Reading graph: %s", graph_path);

  if(!used) cmd_print_usage("Ignored arguments after last --seq");

  if(dump_seq_hist_n > 1) die("Cannot specify --gap-hist <out> more than once");
  if(dump_frag_hist_n > 1) die("Cannot specify --frag-hist <out> more than once");

  // ctx_thread requires output file
  if(!correct_cmd && !args->out_ctp_path)
    cmd_print_usage("--out <out.ctp> is required");

  //
  // Open graph graph file
  //
  GraphFileReader *gfile = &args->gfile;
  graph_file_open(gfile, graph_path);

  if(!correct_cmd && file_filter_into_ncols(&gfile->fltr) > 1)
    die("Please specify a single colour e.g. %s:0", file_filter_path(&gfile->fltr));

  //
  // Open path files
  //
  size_t path_max_usedcols = 0;
  for(i = 0; i < args->gpfiles.len; i++) {
    // file_filter_update_intocol(&args->pfiles.data[i].fltr, 0);
    if(!correct_cmd && file_filter_into_ncols(&args->gpfiles.data[i].fltr) > 1) {
      die("Please specify a single colour e.g. %s:0",
          file_filter_path(&args->gpfiles.data[i].fltr));
    }
    path_max_usedcols = MAX2(path_max_usedcols,
                             file_filter_into_ncols(&args->gpfiles.data[i].fltr));
  }
  args->path_max_usedcols = path_max_usedcols;

  // Check for compatibility between graph files and path files
  graphs_gpaths_compatible(gfile, 1, args->gpfiles.data, args->gpfiles.len, -1);

  // if no paths loaded, set all max_context values to 1, since >1 kmer only
  // useful if can pickup paths
  if(args->gpfiles.len == 0) {
    for(i = 0; i < inputs->len; i++)
      inputs->data[i].crt_params.max_context = 1;
  }

  // Check frag_len_min < frag_len_max
  for(i = 0; i < inputs->len; i++)
  {
    CorrectAlnInput *t = &inputs->data[i];
    t->files.ptr = t;
    if(t->crt_params.frag_len_min > t->crt_params.frag_len_max) {
      die("--min-ins %u is greater than --max-ins %u",
          t->crt_params.frag_len_min, t->crt_params.frag_len_max);
    }
    correct_aln_input_print(&inputs->data[i]);
    args->max_gap_limit = MAX2(args->max_gap_limit, t->crt_params.frag_len_max);
  }

  futil_create_output(args->dump_seq_sizes);
  futil_create_output(args->dump_frag_sizes);
}
