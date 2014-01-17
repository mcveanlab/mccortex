#include "global.h"

#include "seq_file.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "graph_format.h"
#include "loading_stats.h"
#include "build_graph.h"

static const char usage[] =
"usage: "CMD" build [options] <out.ctx>\n"
"  Build a cortex graph.  \n"
"\n"
"  Options:\n"
"    -m <mem>               Memory to use (e.g. 100G or 12M)\n"
"    -k <kmer>              Kmer size\n"
"    --asyncio <N>          Number of input threads\n"
"    --threads <N>          Number of processing threads\n"
"    --sample <name>        Sample name (required before any seq args)\n"
"    --seq <in.fa|fq|sam>   Load sequence data\n"
"    --seq2 <in1> <in2>     Load paired end sequence data\n"
"    --fq_threshold <qual>  Filter quality scores [default: 0 (off)]\n"
"    --fq_offset <offset>   FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"    --cut_hp <len>         Breaks reads at homopolymers >= <len> [default: off]\n"
"    --remove_pcr           Remove (or keep) PCR duplicate reads [default: keep]\n"
"    --keep_pcr             Don't do PCR duplicate removal\n"
"    --load_graph <in.ctx>  Load samples from a graph file (.ctx)\n"
"\n"
"  PCR duplicate removal works by ignoring read pairs (PE-only) if both reads\n"
"  start at the same k-mer as any previous read. Carried out per sample, not \n"
"  per file.\n"
"\n"
"  --sample <name> is required before sequence input can be loaded.\n"
"  Consecutive sequence options are loaded into the same colour.\n"
"\n"
"  --load_graph argument can have colours specifed e.g. in.ctx:0,6-8 will load\n"
"  samples 0,6,7,8.  Graphs are loaded into new colours.\n"
"\n"
" Example:\n"
"   "CMD" build -k 31 -m 60G --fq_threshold 5 --cut_hp 8 \\\n"
"               --sample bob --seq bob.bam \\\n"
"               --sample dylan --remove_pcr --seq dylan.sam \\\n"
"                              --keep_pcr --fq_offset 33 \\\n"
"                              --seq2 dylan.1.fq.gz dylan.2.fq.gz \\\n"
"               bob.dylan.k31.ctx\n"
"\n"
"  See `"CMD" join` to combine .ctx files\n";

// DEV: check we update all ginfo satisfactorily
// DEV: remove CtxBuildInput struct - not really needed
//      keep BuildGraphTask and GraphFileReader in sep arrays,
//      load one lot then the other

typedef struct
{
  size_t colour;
  SeqLoadingStats *stats;
  char *sample_name;
  boolean empty;
} CtxBuildGraphCol;

typedef struct
{
  // Load graph unless it is NULL, then load seq instead
  BuildGraphTask seq_files; // sequence files to be loaded
  GraphFileReader ctx_file; // ctx_file.hdr.num_cols > 0 if set
  CtxBuildGraphCol *col_data;
} CtxBuildInput;

#define input_is_ctx_file(inpt) ((inpt)->ctx_file.hdr.num_of_cols > 0)
#define input_intocol(inpt) ((inpt)->seq_files.colour)
#define input_ncols(inpt) \
        (input_is_ctx_file(inpt) ? graph_file_usedcols(&(inpt)->ctx_file) : 1)

static void print_input(const CtxBuildInput *input)
{
  const BuildGraphTask *task = &input->seq_files;
  Colour colour = input_intocol(input), ncols = input_ncols(input);

  if(input_is_ctx_file(input))
  {
    status("[sample] %zu-%zu: load graph: %s", colour,
           colour+ncols, graph_file_orig_path(&input->ctx_file));
  }
  else
  {
    char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

    if(task->fq_cutoff > 0) sprintf(fqCutoff, "%u", task->fq_cutoff);
    if(task->fq_offset > 0) sprintf(fqOffset, "%u", task->fq_offset);
    if(task->hp_cutoff > 0) sprintf(hpCutoff, "%u", task->hp_cutoff);

    status("[sample] input: %s%s%s; FASTQ offset: %s, threshold: %s; "
           "cut homopolymers: %s; remove PCR duplicates SE: %s, PE: %s;\n",
           task->file1->path,
           task->file2 ? ", " : "", task->file2 ? task->file2->path : "",
           fqOffset, fqCutoff, hpCutoff,
           task->remove_dups_se ? "yes" : "no",
           task->remove_dups_pe ? "yes" : "no");
  }
}

static CtxBuildInput ctx_input_create(char *load_graph_path,
                                      const char *seq_path1, const char *seq_path2,
                                      size_t colour, CtxBuildGraphCol *col_data,
                                      uint32_t fq_offset, uint32_t fq_cutoff,
                                      uint32_t hp_cutoff,
                                      boolean remove_dups_se,
                                      boolean remove_dups_pe,
                                      size_t kmer_size)
{
  assert(colour >= 0);
  assert((load_graph_path == NULL) != (seq_path1 == NULL));
  assert((seq_path2 == NULL) || (seq_path1 != NULL)); // seq2 => seq1
  assert(fq_offset < 128);
  assert(fq_cutoff < 128);
  assert(hp_cutoff < 256);

  seq_file_t *file1 = NULL, *file2 = NULL;
  const char *arg = (seq_path2 == NULL ? "--seq2" : "--seq");

  if(seq_path1 != NULL && (file1 = seq_open(seq_path1)) == NULL)
    die("Cannot read %s file: %s", arg, seq_path1);
  if(seq_path2 != NULL && (file2 = seq_open(seq_path2)) == NULL)
    die("Cannot read %s file: %s", arg, seq_path2);

  BuildGraphTask seq_files = {.file1 = file1, .file2 = file2, .colour = colour,
                              .fq_offset = (uint8_t)fq_offset,
                              .fq_cutoff = (uint8_t)fq_cutoff,
                              .hp_cutoff = (uint8_t)hp_cutoff,
                              .remove_dups_se = remove_dups_se,
                              .remove_dups_pe = remove_dups_pe,
                              .stats = col_data->stats};

  CtxBuildInput input = {.seq_files = seq_files,
                         .ctx_file = INIT_GRAPH_READER_MACRO,
                         .col_data = col_data};

  col_data->empty = false;

  if(load_graph_path != NULL)
  {
    int ret = graph_file_open(&input.ctx_file, load_graph_path, false);

    if(ret == 0)
      print_usage(usage, "Cannot read input graph file: %s", load_graph_path);
    else if(ret < 0)
      print_usage(usage, "Input graph file isn't valid: %s", load_graph_path);

    if(input.ctx_file.hdr.kmer_size != kmer_size) {
      print_usage(usage, "Input graph kmer_size doesn't match [%u vs %zu]",
                  input.ctx_file.hdr.kmer_size, kmer_size);
    }

    input.ctx_file.fltr.intocol = colour;
  }

  return input;
}

static void ctx_input_dealloc(CtxBuildInput *input)
{
  graph_file_dealloc(&input->ctx_file);
}

static void ctx_build_graph_colour_dealloc(CtxBuildGraphCol *col)
{
  seq_loading_stats_free(col->stats);
}

// inputs must be an already allocated array to put inputs into
// Returns index of next untouched argv
// *num_inputs_ptr is the number of inputs
// *num_seq_cols_ptr is the number of colours with sequenced loaded into them
// *num_total_cols_ptr is the total number of colours in the graph
static int load_args(int argc, char **argv, size_t kmer_size,
                     CtxBuildInput *inputs, size_t *num_inputs_ptr,
                     CtxBuildGraphCol *seq_cols, size_t *num_seq_cols_ptr,
                     size_t *num_total_cols_ptr)
{
  int argi;
  size_t num_inputs = 0, num_seq_cols = 0;
  uint32_t fq_offset = 0, fq_cutoff = 0, hp_cutoff = 0;
  // We don't allow remove_dups_se to ever be true currently...
  boolean remove_dups_se = false, remove_dups_pe = false;

  size_t colour = 0;
  boolean sample_named = false, sample_used = false, seq_loaded = false;

  for(argi = 0; argi < argc; argi++)
  {
    if(strcmp(argv[argi],"--fq_threshold") == 0) {
      if(argi + 1 >= argc)
        print_usage(usage, "--fq_threshold <qual> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_offset) || fq_offset > 128)
        die("Invalid --fq_threshold argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_offset") == 0) {
      if(argi + 1 >= argc)
        print_usage(usage, "--fq_offset <offset> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_offset) || fq_offset > 128)
        die("Invalid --fq_offset argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--cut_hp") == 0) {
      if(argi + 1 >= argc)
        print_usage(usage, "--cut_hp <len> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &hp_cutoff))
        die("Invalid --cut_hp argument: %s", argv[argi+1]);
      if(hp_cutoff > UINT8_MAX)
        die("--cut_hp <hp> cannot be greater than %i", UINT8_MAX);
      argi += 1;
    }
    else if(!strcmp(argv[argi],"--remove_pcr")) { remove_dups_pe = true; }
    else if(!strcmp(argv[argi],"--keep_pcr")) { remove_dups_pe = false; }
    else if(strcmp(argv[argi],"--seq") == 0) {
      if(!sample_named)
        print_usage(usage, "Please use --sample <name> before giving sequence");
      if(argi + 1 >= argc)
        print_usage(usage, "--seq <file> requires an argument");

      // Create new task
      inputs[num_inputs++] = ctx_input_create(NULL, argv[argi+1], NULL,
                                              colour, &seq_cols[num_seq_cols-1],
                                              fq_offset, fq_cutoff, hp_cutoff,
                                              remove_dups_se, remove_dups_pe,
                                              kmer_size);
      //
      argi += 1;
      sample_used = true;
      seq_loaded = true;
    }
    else if(strcmp(argv[argi],"--seq2") == 0) {
      if(!sample_named)
        print_usage(usage, "Please use --sample <name> before giving sequence");
      if(argi + 2 >= argc)
        print_usage(usage, "--seq2 <file1> <file2> requires two arguments");

      // Create new task
      inputs[num_inputs++] = ctx_input_create(NULL, argv[argi+1], argv[argi+2],
                                              colour, &seq_cols[num_seq_cols-1],
                                              fq_offset, fq_cutoff, hp_cutoff,
                                              remove_dups_se, remove_dups_pe,
                                              kmer_size);
      //
      argi += 2;
      sample_used = true;
      seq_loaded = true;
    }
    else if(!strcmp(argv[argi],"--sample") || !strcmp(argv[argi],"--load_graph"))
    {
      if(argi + 1 >= argc) print_usage(usage, "%s requires an arg", argv[argi]);

      if(sample_named && !sample_used) {
        warn("Empty colour '%s' (maybe you intended this?)",
             seq_cols[num_seq_cols-1].sample_name);
      }

      if(sample_named) { colour++; sample_named = sample_used = false; }

      if(!strcmp(argv[argi],"--sample"))
      {
        // --sample <name>
        if(!argv[argi+1][0] || !strcmp(argv[argi+1], "undefined") ||
           !strcmp(argv[argi+1],"noname")) {
          die("--sample %s is not a good name!", argv[argi+1]);
        }

        seq_cols[num_seq_cols].colour = colour;
        seq_cols[num_seq_cols].stats = seq_loading_stats_create(0);
        seq_cols[num_seq_cols].sample_name = argv[argi+1];
        seq_cols[num_seq_cols].empty = true;
        num_seq_cols++;

        sample_named = true;
      }
      else
      {
        // --load_graph <in.ctx>
        // Load binary into new colour
        inputs[num_inputs] = ctx_input_create(argv[argi+1], NULL, NULL,
                                              colour, NULL,
                                              fq_offset, fq_cutoff, hp_cutoff,
                                              remove_dups_se, remove_dups_pe,
                                              kmer_size);
        colour += input_ncols(&inputs[num_inputs]);
        num_inputs++;
      }

      argi++; // both --sample and --load_graph take a single argument
    }
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  if(!seq_loaded) die("No sequence loaded, to combine graph files use '"CMD" join'");
  if(sample_named && !sample_used) {
    warn("Empty colour '%s' (maybe you intended this?)",
         seq_cols[num_seq_cols-1].sample_name);
  }
  if(sample_named) colour++;

  *num_inputs_ptr = num_inputs;
  *num_seq_cols_ptr = num_seq_cols;
  *num_total_cols_ptr = colour;

  return argi;
}

static void run_build_graph(dBGraph *db_graph,
                            const CtxBuildInput *inputs, size_t num_inputs,
                            size_t num_build_threads)
{
  status("Loading %zu files with %zu worker threads", num_inputs, num_build_threads);
  size_t i;
  BuildGraphTask *tasks = malloc(num_inputs*sizeof(BuildGraphTask));
  for(i = 0; i < num_inputs; i++) tasks[i] = inputs[i].seq_files;
  build_graph(db_graph, tasks, num_inputs, num_build_threads);
  free(tasks);
}

int ctx_build(CmdArgs *args)
{
  cmd_accept_options(args, "atmnk", usage);
  cmd_require_options(args, "k", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  size_t kmer_size = args->kmer_size;

  // num_seq_cols is the number of colours with sequence loaded into them
  // output_colours = num_seq_cols + num colours filled with graph files
  size_t i, num_inputs = 0, num_seq_cols = 0, output_colours = 0;
  int argi = 0;

  // Load inputs
  size_t max_inputs = (size_t)argc/2;
  CtxBuildInput *inputs = malloc2(max_inputs * sizeof(CtxBuildInput));
  CtxBuildGraphCol *seq_cols = malloc2(max_inputs * sizeof(CtxBuildGraphCol));

  argi = load_args(argc-1, argv, kmer_size, inputs, &num_inputs,
                   seq_cols, &num_seq_cols, &output_colours);

  const char *out_path = argv[argc-1];

  // Did any inputs require PCR duplicate removal
  boolean remove_pcr_used = false;

  for(i = 0; i < num_inputs; i++) {
    if(input_is_ctx_file(&inputs[i]) &&
       (inputs[i].seq_files.remove_dups_se || inputs[i].seq_files.remove_dups_pe))
    {
      remove_pcr_used = true;
      break;
    }
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = ((sizeof(Covg) + sizeof(Edges)) * 8 + remove_pcr_used*2) *
                  output_colours;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, 0, true, &graph_mem);
  cmd_check_mem_limit(args, graph_mem);

  //
  // Check output path
  //
  status("Writing %zu colour graph to %s\n", output_colours, out_path);

  if(!futil_is_file_writable(out_path))
    die("Cannot write to file: %s", out_path);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, output_colours, output_colours, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity * output_colours, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity * output_colours, sizeof(Covg));

  size_t bktlocks_mem = roundup_bits2bytes(db_graph.ht.num_of_buckets);
  db_graph.bktlocks = calloc2(bktlocks_mem, sizeof(uint8_t));

  size_t kmer_words = roundup_bits2words64(db_graph.ht.capacity);

  if(remove_pcr_used) {
    db_graph.readstrt = calloc2(roundup_bits2bytes(db_graph.ht.capacity)*2,
                                sizeof(uint8_t));
  }

  hash_table_print_stats(&db_graph.ht);

  size_t j, colour, ncols, from_col;
  GraphFileReader *ctx_file;

  // Set sample names using seq_colours array
  for(i = 0; i < num_seq_cols; i++) {
    strbuf_set(&db_graph.ginfo[seq_cols[i].colour].sample_name,
               seq_cols[i].sample_name);
  }

  for(i = 0; i < num_inputs; i++) {
    if(i == 0 || colour != input_intocol(&inputs[i])) {
      // New colour
      colour = input_intocol(&inputs[i]);

      if(input_is_ctx_file(&inputs[i])) {
        ncols = input_ncols(&inputs[i]);
        ctx_file = &inputs[i].ctx_file;
        status("[sample] %zu-%zu: %s\n", colour, colour + ncols - 1,
               graph_file_orig_path(ctx_file));
        for(j = 0; j < ncols; j++) {
          from_col = graph_file_fromcol(ctx_file, j);
          status("[sample]   %s\n", ctx_file->hdr.ginfo[from_col].sample_name.buff);
        }
      }
      else {
        status("[sample] %zu: %s\n", colour, inputs[i].col_data->sample_name);
      }
    }
    print_input(&inputs[i]);
  }

  size_t start, end, num_load, prev_colour = 0;

  // If we are using PCR duplicate removal,
  // best to load one colour at a time
  for(start = 0; start < num_inputs; start = end, prev_colour = colour)
  {
    // Wipe read starts
    colour = input_intocol(&inputs[start]);
    if(remove_pcr_used)
    {
      if(colour != prev_colour)
        memset(db_graph.readstrt, 0, 2*kmer_words*sizeof(uint64_t));

      end = start+1;
      while(end < num_inputs && end-start < args->max_io_threads &&
            input_intocol(&inputs[end]) == colour) end++;
    }
    else end = MIN2(start+args->max_io_threads, num_inputs);

    num_load = end-start;
    run_build_graph(&db_graph, inputs+start, num_load, args->max_work_threads);
  }

  hash_table_print_stats(&db_graph.ht);

  // Update ginfo
  for(i = 0; i < num_seq_cols; i++) {
    colour = seq_cols[i].colour;
    graph_info_update_contigs(&db_graph.ginfo[colour],
                              seq_cols[i].stats->total_bases_loaded,
                              seq_cols[i].stats->contigs_loaded);
  }

  status("Dumping graph...\n");
  graph_file_save_mkhdr(out_path, &db_graph, CTX_GRAPH_FILEFORMAT, NULL,
                        0, output_colours);

  for(i = 0; i < num_inputs; i++) ctx_input_dealloc(&inputs[i]);
  for(i = 0; i < num_seq_cols; i++) ctx_build_graph_colour_dealloc(&seq_cols[i]);

  free(inputs);
  free(seq_cols);

  free((uint8_t*)db_graph.bktlocks);
  free(db_graph.col_covgs);
  free(db_graph.col_edges);
  if(db_graph.readstrt != NULL) free(db_graph.readstrt);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
