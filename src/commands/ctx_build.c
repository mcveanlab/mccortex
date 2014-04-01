#include "global.h"

#include "seq_file.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "graph_format.h"
#include "loading_stats.h"
#include "build_graph.h"

const char build_usage[] =
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
"    --seqi <in.bam>        Load paired end sequence from a single file\n"
"    --fq_threshold <qual>  Filter quality scores [default: 0 (off)]\n"
"    --fq_offset <offset>   FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"    --cut_hp <len>         Breaks reads at homopolymers >= <len> [default: off]\n"
"    --remove_pcr           Remove (or keep) PCR duplicate reads [default: keep]\n"
"    --keep_pcr             Don't do PCR duplicate removal\n"
"    --load_graph <in.ctx>  Load samples from a graph file (.ctx)\n"
"    --FR --FF --RF --RR    Mate pair orientation [default: FR] (used with --keep_pcr)\n"
"\n"
"  PCR duplicate removal works by ignoring read (pairs) if (both) reads\n"
"  start at the same k-mer as any previous read. Carried out per sample, not \n"
"  per file. --sample <name> is required before sequence input can be loaded.\n"
"  Consecutive sequence options are loaded into the same colour.\n"
"  --load_graph argument can have colours specifed e.g. in.ctx:0,6-8 will load\n"
"  samples 0,6,7,8.  Graphs are loaded into new colours.\n"
"  See `"CMD" join` to combine .ctx files\n";

typedef struct
{
  size_t colour;
  char *sample_name;
} SampleName;

static void print_graph_reader(const GraphFileReader *rdr)
{
  Colour colour = graph_file_intocol(rdr,0), ncols = graph_file_outncols(rdr);
  const char *path = graph_file_orig_path(rdr);

  status("[load] %zu-%zu: load graph: %s", colour, colour+ncols-1, path);
}

static void graph_reader_new(char *path, size_t colour, GraphFileReader *rdrptr)
{
  GraphFileReader rdr = INIT_GRAPH_READER_MACRO;
  int ret = graph_file_open(&rdr, path, false);

  if(ret == 0) cmd_print_usage("Cannot read input graph file: %s", path);
  else if(ret < 0) cmd_print_usage("Input graph file isn't valid: %s", path);

  rdr.fltr.intocol = colour;

  memcpy(rdrptr, &rdrptr, sizeof(GraphFileReader));
}

static void build_graph_task_new(const char *seq_path1,
                                 const char *seq_path2,
                                 bool interleaved,
                                 size_t colour,
                                 uint32_t fq_offset,
                                 uint32_t fq_cutoff,
                                 uint32_t hp_cutoff,
                                 bool remove_pcr_dups,
                                 ReadMateDir matedir,
                                 BuildGraphTask *taskptr)
{
  ctx_assert(!(interleaved && seq_path2));
  ctx_assert((seq_path2 == NULL) || (seq_path1 != NULL)); // seq2 => seq1
  ctx_assert(fq_offset < 128);
  ctx_assert(fq_cutoff < 128);
  ctx_assert(hp_cutoff < 256);

  seq_file_t *file1 = NULL, *file2 = NULL;
  const char *arg = (seq_path2 == NULL ? "--seq" : "--seq2");

  if(seq_path1 != NULL && (file1 = seq_open(seq_path1)) == NULL)
    die("Cannot read %s file: %s", arg, seq_path1);
  if(seq_path2 != NULL && (file2 = seq_open(seq_path2)) == NULL)
    die("Cannot read %s file: %s", arg, seq_path2);

  AsyncIOReadTask iotask = {.file1 = file1, .file2 = file2,
                            .fq_offset = (uint8_t)fq_offset,
                            .interleaved = interleaved,
                            .ptr = taskptr};

  BuildGraphTask task = {.files = iotask,
                         .fq_cutoff = (uint8_t)fq_cutoff,
                         .hp_cutoff = (uint8_t)hp_cutoff,
                         .remove_pcr_dups = remove_pcr_dups,
                         .matedir = matedir,
                         .colour = colour,
                         .stats = LOAD_STATS_INIT_MACRO};

  memcpy(taskptr, &task, sizeof(BuildGraphTask));
}

static void build_graph_task_destroy(BuildGraphTask *task)
{
  asyncio_task_close(&task->files);
}

static void build_graph_task_new2(const char *seq_path1, const char *seq_path2,
                                  bool interleaved, size_t colour,
                                  uint32_t fq_offset, uint32_t fq_cutoff,
                                  uint32_t hp_cutoff, bool remove_pcr_dups,
                                  ReadMateDir matedir,
                                  BuildGraphTask *tasks, size_t *num_tasks_ptr)
{
  ctx_assert(!seq_path2 || seq_path1);
  ctx_assert(!(interleaved && seq_path2));

  if(remove_pcr_dups) {
    // Submit paired end reads together
    build_graph_task_new(seq_path1, seq_path2, interleaved, colour,
                         fq_offset, fq_cutoff, hp_cutoff,
                         remove_pcr_dups, matedir,
                         &tasks[*num_tasks_ptr]);
    (*num_tasks_ptr)++;
  }
  else if(!remove_pcr_dups) {
    // Read files separately -> read faster
    build_graph_task_new(seq_path1, NULL, interleaved, colour,
                         fq_offset, fq_cutoff, hp_cutoff,
                         remove_pcr_dups, matedir,
                         &tasks[*num_tasks_ptr]);
    (*num_tasks_ptr)++;

    if(seq_path2 != NULL) {
      build_graph_task_new(seq_path2, NULL, false, colour,
                           fq_offset, fq_cutoff, hp_cutoff,
                           remove_pcr_dups, matedir,
                           &tasks[*num_tasks_ptr]);
      (*num_tasks_ptr)++;
    }
  }
}

static void load_args(int argc, char **argv, size_t *kmer_size_ptr,
                      BuildGraphTask *tasks, size_t *num_tasks_ptr,
                      GraphFileReader *graphs, size_t *num_graphs_ptr,
                      SampleName *samples, size_t *num_samples_ptr,
                      size_t *num_total_cols_ptr)
{
  int argi;
  size_t num_graphs = 0, num_tasks = 0, num_samples = 0;
  uint32_t fq_offset = 0, fq_cutoff = 0, hp_cutoff = 0;
  bool remove_pcr_dups = false;
  ReadMateDir matedir = READPAIR_FR;

  size_t i, colour = SIZE_MAX, kmer_size = 0;
  bool sample_named = false, sample_used = false, seq_loaded = false;

  for(argi = 0; argi < argc; argi++)
  {
    if(!strcmp(argv[argi],"--kmer_size") || !strcmp(argv[argi],"-k")) {
      if(argi + 1 >= argc) die("%s <k> requires an argument", argv[argi]);
      if(!parse_entire_size(argv[argi+1], &kmer_size) || !(kmer_size&1)) {
        die("Invalid kmer-size (%s %s): requires odd number %i <= k <= %i",
            argv[argi], argv[argi+1], MIN_KMER_SIZE, MAX_KMER_SIZE);
      }
      if(kmer_size < MIN_KMER_SIZE || kmer_size > MAX_KMER_SIZE) {
        die("Please recompile with correct kmer size (%zu)", kmer_size);
      }
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_threshold") == 0) {
      if(argi + 1 >= argc)
        cmd_print_usage("--fq_threshold <qual> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_cutoff) || fq_cutoff > 128)
        die("Invalid --fq_threshold argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_offset") == 0) {
      if(argi + 1 >= argc)
        cmd_print_usage("--fq_offset <offset> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_offset) || fq_offset > 128)
        die("Invalid --fq_offset argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--cut_hp") == 0) {
      if(argi + 1 >= argc)
        cmd_print_usage("--cut_hp <len> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &hp_cutoff))
        die("Invalid --cut_hp argument: %s", argv[argi+1]);
      if(hp_cutoff > UINT8_MAX)
        die("--cut_hp <hp> cannot be greater than %i", UINT8_MAX);
      argi += 1;
    }
    else if(strcasecmp(argv[argi],"--FF") == 0) matedir = READPAIR_FF;
    else if(strcasecmp(argv[argi],"--FR") == 0) matedir = READPAIR_FR;
    else if(strcasecmp(argv[argi],"--RF") == 0) matedir = READPAIR_RF;
    else if(strcasecmp(argv[argi],"--RR") == 0) matedir = READPAIR_RR;
    else if(!strcmp(argv[argi],"--remove_pcr")) remove_pcr_dups = true;
    else if(!strcmp(argv[argi],"--keep_pcr")) remove_pcr_dups = false;
    else if(strcmp(argv[argi],"--seq") == 0 || strcmp(argv[argi],"--seqi") == 0)
    {
      if(!sample_named)
        cmd_print_usage("Please use --sample <name> before giving sequence");
      if(argi + 1 >= argc)
        cmd_print_usage("--seq <file> requires an argument");

      bool interleaved = (strcmp(argv[argi],"--seqi") == 0);

      // Create new task
      build_graph_task_new2(argv[argi+1], NULL, interleaved, colour,
                            fq_offset, fq_cutoff, hp_cutoff,
                            remove_pcr_dups, matedir,
                            tasks, &num_tasks);
      //
      argi += 1;
      sample_used = true;
      seq_loaded = true;
    }
    else if(strcmp(argv[argi],"--seq2") == 0) {
      if(!sample_named)
        cmd_print_usage("Please use --sample <name> before giving sequence");
      if(argi + 2 >= argc)
        cmd_print_usage("--seq2 <file1> <file2> requires two arguments");

      // Create new task
      build_graph_task_new2(argv[argi+1], argv[argi+2], false, colour,
                            fq_offset, fq_cutoff, hp_cutoff,
                            remove_pcr_dups, matedir,
                            tasks, &num_tasks);
      //
      argi += 2;
      sample_used = true;
      seq_loaded = true;
    }
    else if(!strcmp(argv[argi],"--sample") || !strcmp(argv[argi],"--load_graph"))
    {
      if(argi + 1 >= argc) cmd_print_usage("%s requires an arg", argv[argi]);

      if(sample_named && !sample_used) {
        warn("Empty colour '%s' (maybe you intended this?)",
             samples[num_samples-1].sample_name);
      }

      // Move on to next colour
      colour++;

      if(strcmp(argv[argi],"--sample") == 0)
      {
        // --sample <name>
        if(!argv[argi+1][0] || !strcmp(argv[argi+1], "undefined") ||
           !strcmp(argv[argi+1],"noname")) {
          die("--sample %s is not a good name!", argv[argi+1]);
        }

        samples[num_samples].colour = colour;
        samples[num_samples].sample_name = argv[argi+1];
        num_samples++;

        sample_named = true;
        sample_used = false;
      }
      else
      {
        // --load_graph <in.ctx>
        // Load binary into new colour
        graph_reader_new(argv[argi+1], colour, &graphs[num_graphs]);
        colour += graph_file_outncols(&graphs[num_graphs]);
        num_graphs++;
        sample_named = false;
        sample_used = false;
      }

      argi++; // both --sample and --load_graph take a single argument
    }
    else cmd_print_usage("Unknown option: %s", argv[argi]);
  }

  if(!seq_loaded)
    die("No sequence loaded, to combine graph files use '"CMD" join'");

  if(sample_named && !sample_used) {
    warn("Empty colour '%s' (maybe you intended this?)",
         samples[num_samples-1].sample_name);
  }

  if(kmer_size == 0)
    die("--kmer_size <k> not set");

  if(sample_named) colour++;

  // Check kmer size in graphs to load
  for(i = 0; i < num_graphs; i++) {
    if(graphs[i].hdr.kmer_size != kmer_size) {
      cmd_print_usage("Input graph kmer_size doesn't match [%u vs %zu]: %s",
                      graphs[i].hdr.kmer_size, kmer_size,
                      graphs[i].fltr.orig_path.buff);
    }
  }

  *num_graphs_ptr = num_graphs;
  *num_tasks_ptr = num_tasks;
  *num_samples_ptr = num_samples;
  *num_total_cols_ptr = colour;
  *kmer_size_ptr = kmer_size;
}

int ctx_build(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that we have at least 1 argument

  // num_seq_cols is the number of colours with sequence loaded into them
  // output_colours = num_seq_cols + num colours filled with graph files
  size_t i, kmer_size;
  size_t num_tasks = 0, num_graphs = 0, num_samples = 0, output_colours = 0;

  // Load inputs
  size_t max_inputs = (size_t)argc/2;

  BuildGraphTask *tasks = malloc2(max_inputs * sizeof(BuildGraphTask));
  GraphFileReader *graphs = malloc2(max_inputs * sizeof(GraphFileReader));
  SampleName *samples = malloc2(max_inputs * sizeof(GraphFileReader));

  load_args(argc-1, argv, &kmer_size,
            tasks, &num_tasks,
            graphs, &num_graphs,
            samples, &num_samples,
            &output_colours);

  const char *out_path = argv[argc-1];

  // Did any tasks require PCR duplicate removal
  for(i = 0; i < num_tasks && !tasks[i].remove_pcr_dups; i++) {}
  bool remove_pcr_used = (i < num_tasks);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  // remove_pcr_dups requires a fw and rv bit per kmer
  bits_per_kmer = (sizeof(Covg) + sizeof(Edges))*8*output_colours +
                  remove_pcr_used*2;

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, 0, true, &graph_mem);
  cmd_check_mem_limit(args, graph_mem);

  //
  // Check output path
  //
  const char *out_path_name = (strcmp(out_path,"-") == 0 ? "STDOUT" : out_path);

  if(strcmp(out_path,"-") != 0)
  {
    if(futil_file_exists(out_path)) die("Output file already exists: %s", out_path);
    if(!futil_is_file_writable(out_path)) die("Cannot write to file: %s", out_path);
  }

  // Print graphs to be loaded
  for(i = 0; i < num_graphs; i++)
    print_graph_reader(&graphs[i]);

  // Print tasks and sample names
  size_t s;
  status("[sample] %zu: %s", (size_t)0, samples[0].sample_name);

  for(i = 0, s = 0; i < num_tasks; i++) {
    while(samples[s].colour < tasks[i].colour) {
      s++; status("[sample] %zu: %s", s, samples[s].sample_name);
    }
    build_graph_task_print(&tasks[i]);
  }

  // Print remaining empty samples
  for(; s < num_samples; s++) {
    status("[sample] %zu: %s", s, samples[s].sample_name);
  }

  status("Writing %zu colour graph to %s\n", output_colours, out_path_name);

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

  // Load graphs
  if(num_graphs > 0)
  {
    GraphLoadingPrefs gprefs = LOAD_GPREFS_INIT(&db_graph);
    LoadingStats gstats = LOAD_STATS_INIT_MACRO;

    for(i = 0; i < num_graphs; i++) {
      graph_load(&graphs[i], gprefs, &gstats);
      graph_file_dealloc(&graphs[i]);
    }

    hash_table_print_stats(&db_graph.ht);
  }

  // Set sample names using seq_colours array
  for(i = 0; i < num_samples; i++) {
    strbuf_set(&db_graph.ginfo[samples[i].colour].sample_name,
               samples[i].sample_name);
  }

  size_t start, end, num_load, colour, prev_colour = 0;

  // If we are using PCR duplicate removal,
  // it's best to load one colour at a time
  for(start = 0; start < num_tasks; start = end, prev_colour = colour)
  {
    // Wipe read starts
    colour = tasks[start].colour;
    if(remove_pcr_used)
    {
      if(colour != prev_colour)
        memset(db_graph.readstrt, 0, 2*kmer_words*sizeof(uint64_t));

      end = start+1;
      while(end < num_tasks && end-start < args->max_io_threads &&
            tasks[end].colour == colour) end++;
    }
    else {
      end = MIN2(start+args->max_io_threads, num_tasks);
    }

    num_load = end-start;
    build_graph(&db_graph, tasks+start, num_load, args->max_work_threads);
  }

  // Print stats for hash table
  hash_table_print_stats(&db_graph.ht);

  // Print stats per input file
  for(i = 0; i < num_tasks; i++) {
    build_graph_task_print_stats(&tasks[i]);
    build_graph_task_destroy(&tasks[i]);
  }

  status("Dumping graph...\n");
  graph_file_save_mkhdr(out_path, &db_graph, CTX_GRAPH_FILEFORMAT, NULL,
                        0, output_colours);

  free(tasks);
  free(samples);
  free(graphs);

  free(db_graph.bktlocks);
  free(db_graph.col_covgs);
  free(db_graph.col_edges);
  if(db_graph.readstrt != NULL) free(db_graph.readstrt);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
