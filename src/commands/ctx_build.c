#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "graphs_load.h"
#include "graph_writer.h"
#include "build_graph.h"

#include "seq_file/seq_file.h"

const char build_usage[] =
"usage: "CMD" build [options] <out.ctx>\n"
"\n"
"  Build a cortex graph.  \n"
"\n"
"  -h, --help               This help message\n"
"  -q, --quiet              Silence status output normally printed to STDERR\n"
"  -f, --force              Overwrite output files\n"
"  -m, --memory <mem>       Memory to use\n"
"  -n, --nkmers <kmers>     Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>        Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
//
"  -k, --kmer <kmer>        Kmer size must be odd ("QUOTE_VALUE(MAX_KMER_SIZE)" >= k >= "QUOTE_VALUE(MIN_KMER_SIZE)")\n"
"  -s, --sample <name>      Sample name (required before any seq args)\n"
"  -1, --seq <in.fa>        Load sequence data\n"
"  -2, --seq2 <in1:in2>     Load paired end sequence data\n"
"  -i, --seqi <in.bam>      Load paired end sequence from a single file\n"
"  -Q, --fq-cutoff <Q>      Filter quality scores [default: 0 (off)]\n"
"  -O, --fq-offset <N>      FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"  -H, --cut-hp <bp>        Breaks reads at homopolymers >= <bp> [default: off]\n"
"  -p, --remove-pcr         Remove (or keep) PCR duplicate reads\n"
"  -P, --keep-pcr           Don't do PCR duplicate removal [default]\n"
"  -M, --matepair <orient>  Mate pair orientation: FF,FR,RF,RR [default: FR]\n"
"                           (for --keep_pcr only)\n"
"  -g, --graph <in.ctx>     Load samples from a graph file (.ctx)\n"
"  -I, --intersect <i.ctx>  Only load kmers that appear in i.ctx. Multiple -I\n"
"                           graphs will be merged, not intersected. Treated as\n"
"                           single colour graphs.\n"
" -S, --sort                Output a graph file ordered by kmer\n"
"\n"
"  Note: Argument must come before input file\n"
"  PCR duplicate removal works by ignoring read (pairs) if (both) reads\n"
"  start at the same k-mer as any previous read. Carried out per sample, not \n"
"  per file. --sample <name> is required before sequence input can be loaded.\n"
"  Consecutive sequence options are loaded into the same colour.\n"
"  --graph argument can have colours specifed e.g. in.ctx:0,6-8 will load\n"
"  samples 0,6,7,8.  Graphs are loaded into new colours.\n"
"  See `"CMD" join` to combine .ctx files\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"kmer",         required_argument, NULL, 'k'},
  {"sample",       required_argument, NULL, 's'},
  {"sort",         no_argument,       NULL, 'S'},
  {"seq",          required_argument, NULL, '1'},
  {"seq2",         required_argument, NULL, '2'},
  {"seqi",         required_argument, NULL, 'i'},
  {"matepair",     required_argument, NULL, 'M'},
  {"fq-cutoff",    required_argument, NULL, 'Q'},
  {"fq-offset",    required_argument, NULL, 'O'},
  {"cut-hp",       required_argument, NULL, 'H'},
  {"remove-pcr",   no_argument,       NULL, 'p'},
  {"keep-pcr",     no_argument,       NULL, 'P'},
  {"graph",        required_argument, NULL, 'g'},
  {"intersect",    required_argument, NULL, 'I'},
  {NULL, 0, NULL, 0}
};

typedef struct {
  size_t colour;
  const char *name;
} SampleName;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(sample_name_buf, SampleNameBuffer, SampleName);

static BuildGraphTaskBuffer gtaskbuf;
static GraphFileBuffer gfilebuf, gisecbuf;
static SampleNameBuffer snamebuf;

static size_t nthreads = 0;
static struct MemArgs memargs = MEM_ARGS_INIT;

static char *out_path = NULL;
static size_t output_colours = 0, kmer_size = 0;

static bool sort_kmers = false;

static void add_task(BuildGraphTask *task)
{
  uint8_t fq_offset = task->files.fq_offset, fq_cutoff = task->prefs.fq_cutoff;
  if(fq_offset >= 128) die("fq-offset too big: %i", (int)fq_offset);
  if(fq_offset+fq_cutoff >= 128) die("fq-cutoff too big: %i", fq_offset+fq_cutoff);

  if(task->prefs.remove_pcr_dups || task->files.file2 == NULL) {
    // Submit paired end reads together
    build_graph_task_buf_push(&gtaskbuf, task, 1);
  }
  else {
    // Read files separately -> read faster
    BuildGraphTask task2 = *task;
    task2.files.file1 = task->files.file2;
    task->files.file2 = task2.files.file2 = NULL;
    build_graph_task_buf_push(&gtaskbuf, task, 1);
    build_graph_task_buf_push(&gtaskbuf, &task2, 1);
  }
}

// Exit with error if bad sample name
static void check_sample_name(const char *sname)
{
  const char *s;
  if(strlen(sname) < 1) die("Sample name is too short: '%s'", sname);
  if(strcmp(sname,"undefined") == 0) die("Bad sample name: '%s'", sname);
  if(strcmp(sname,"noname") == 0) die("Bad sample name: '%s'", sname);
  if(sname[0] == '.') die("Sample name should start with a dot: '%s'", sname);
  for(s = sname; *s; s++) {
    if(isspace(*s)) die("Sample name should not contain whitespace: '%s'", sname);
    if(!isgraph(*s)) die("Bad character in sample name: '%s'", sname);
  }
}

static void parse_args(int argc, char **argv)
{
  BuildGraphTask task;
  memset(&task, 0, sizeof(task));
  task.prefs = SEQ_LOADING_PREFS_INIT;
  task.stats = SEQ_LOADING_STATS_INIT;
  uint8_t fq_offset = 0;
  int intocolour = -1;
  GraphFileReader tmp_gfile;

  // Arg parsing
  char cmd[100], shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;
  bool sample_named = false, pref_unused = false;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 't': cmd_check(!nthreads,cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'k': cmd_check(!kmer_size,cmd); kmer_size = cmd_kmer_size(cmd, optarg); break;
      case 's':
        intocolour++;
        check_sample_name(optarg);
        sample_name_buf_add(&snamebuf, (SampleName){.colour = intocolour,
                                                    .name = optarg});
        sample_named = true;
        break;
      case 'S': cmd_check(!sort_kmers,cmd); sort_kmers = true; break;
      case '1':
      case '2':
      case 'i':
        pref_unused = false;
        if(!sample_named)
          cmd_print_usage("Please give sample name first [-s,--sample <name>]");
        asyncio_task_parse(&task.files, c, optarg, fq_offset, NULL);
        task.prefs.colour = intocolour;
        add_task(&task);
        break;
      case 'M':
             if(!strcmp(optarg,"FF")) task.prefs.matedir = READPAIR_FF;
        else if(!strcmp(optarg,"FR")) task.prefs.matedir = READPAIR_FR;
        else if(!strcmp(optarg,"RF")) task.prefs.matedir = READPAIR_RF;
        else if(!strcmp(optarg,"RR")) task.prefs.matedir = READPAIR_RR;
        else die("-M,--matepair <orient> must be one of: FF,FR,RF,RR");
        pref_unused = true; break;
      case 'O': fq_offset = cmd_uint8(cmd, optarg); pref_unused = true; break;
      case 'Q': task.prefs.fq_cutoff = cmd_uint8(cmd, optarg); pref_unused = true; break;
      case 'H': task.prefs.hp_cutoff = cmd_uint8(cmd, optarg); pref_unused = true; break;
      case 'p': task.prefs.remove_pcr_dups = true; pref_unused = true; break;
      case 'P': task.prefs.remove_pcr_dups = false; pref_unused = true; break;
      case 'g':
        if(intocolour == -1) intocolour = 0;
        graph_file_reset(&tmp_gfile);
        graph_file_open2(&tmp_gfile, optarg, "r", true, intocolour);
        intocolour = MAX2((size_t)intocolour, file_filter_into_ncols(&tmp_gfile.fltr)-1);
        gfile_buf_push(&gfilebuf, &tmp_gfile, 1);
        sample_named = false;
        break;
      case 'I':
        graph_file_reset(&tmp_gfile);
        graph_file_open(&tmp_gfile, optarg);
        if(file_filter_into_ncols(&tmp_gfile.fltr) > 1)
          warn("Flattening intersection graph into colour 0: %s", optarg);
        file_filter_flatten(&tmp_gfile.fltr, 0);
        gfile_buf_push(&gisecbuf, &tmp_gfile, 1);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" build -h` for help. Bad option: %s", argv[optind-1]);
      default: die("Bad option: %s", cmd);
    }
  }

  // Defaults
  if(!nthreads) nthreads = DEFAULT_NTHREADS;

  // Check that optind+1 == argc
  if(optind+1 > argc)
    cmd_print_usage("Expected exactly one graph file");
  else if(optind+1 < argc)
    cmd_print_usage("Expected only one graph file. What is this: '%s'", argv[optind]);

  out_path = argv[optind];
  status("Saving graph to: %s", futil_outpath_str(out_path));

  if(snamebuf.len == 0) cmd_print_usage("No inputs given");

  if(pref_unused) cmd_print_usage("Arguments not given BEFORE sequence file");

  if(!kmer_size) die("kmer size not set with -k <K>");

  // Check kmer size in graphs to load
  size_t i;
  for(i = 0; i < gfilebuf.len; i++) {
    if(gfilebuf.b[i].hdr.kmer_size != kmer_size) {
      cmd_print_usage("Input graph kmer_size doesn't match [%u vs %zu]: %s",
                      gfilebuf.b[i].hdr.kmer_size, kmer_size,
                      file_filter_input(&gfilebuf.b[i].fltr));
    }
  }

  output_colours = intocolour + (sample_named ? 1 : 0);
}


int ctx_build(int argc, char **argv)
{
  size_t i;
  build_graph_task_buf_alloc(&gtaskbuf, 16);
  gfile_buf_alloc(&gfilebuf, 8);
  gfile_buf_alloc(&gisecbuf, 8);
  sample_name_buf_alloc(&snamebuf, 16);

  parse_args(argc, argv);

  size_t s, t, ncolours = snamebuf.len, ntasks = gtaskbuf.len;
  SampleName *samples = snamebuf.b;
  BuildGraphTask *tasks = gtaskbuf.b;

  // Did any tasks require PCR duplicate removal
  for(i = 0; i < ntasks && !tasks[i].prefs.remove_pcr_dups; i++) {}
  bool remove_pcr_used = (i < ntasks);

  //
  // Print inputs
  //
  size_t max_kmers = 0;

  // Print graphs to be loaded
  for(i = 0; i < gfilebuf.len; i++) {
    file_filter_status(&gfilebuf.b[i].fltr);
    max_kmers += gfilebuf.b[i].num_of_kmers;
  }

  // Print tasks and sample names
  for(s = t = 0; s < ncolours || t < ntasks; ) {
    if(t == ntasks || (s < ncolours && samples[s].colour <= tasks[t].prefs.colour)) {
      status("[sample] %zu: %s", s, samples[s].name);
      s++;
    } else {
      build_graph_task_print(&tasks[t]);
      t++;
    }
  }

  for(t = 0; t < ntasks; t++) {
    size_t nkmers = asyncio_input_nkmers(&tasks[t].files);
    if(nkmers == SIZE_MAX) { max_kmers = nkmers; break; }
    max_kmers += nkmers;
  }

  // Check if we are intersecting with graphs
  if(gisecbuf.len > 0)
  {
    if(remove_pcr_used)
      cmd_print_usage("Cannot use --remove-pcr and --intersect");

    for(t = 0; t < ntasks; t++)
      tasks[t].prefs.must_exist_in_graph = true;

    max_kmers = 0;
    for(i = 0; i < gisecbuf.len; i++)
      max_kmers += gisecbuf.b[i].num_of_kmers;
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  // remove_pcr_dups requires a fw and rv bit per kmer
  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  (sizeof(Covg) + sizeof(Edges)) * 8 * output_colours +
                  (gisecbuf.len > 0 ? sizeof(Edges)*8 : 0) +
                  (remove_pcr_used ? 2 : 0) +
                  (sort_kmers ? sizeof(hkey_t) : 0);

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer, 0, max_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Check output path
  //
  futil_create_output(out_path);

  status("Writing %zu colour graph to %s\n", output_colours, futil_outpath_str(out_path));

  // Create db_graph
  dBGraph db_graph;
  int alloc_flags = DBG_ALLOC_EDGES | DBG_ALLOC_COVGS | DBG_ALLOC_BKTLOCKS |
                    (remove_pcr_used ? DBG_ALLOC_READSTRT : 0);

  db_graph_alloc(&db_graph, kmer_size, output_colours, output_colours,
                 kmers_in_hash, alloc_flags);

  Edges *isec_edges = NULL;
  if(gisecbuf.len > 0)
    isec_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));

  hash_table_print_stats(&db_graph.ht);

  // Load intersection graphs
  if(gisecbuf.len > 0)
  {
    GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);
    Covg *tmp_covgs = NULL;
    SWAP(db_graph.col_covgs, tmp_covgs);
    SWAP(db_graph.col_edges, isec_edges); db_graph.num_edge_cols = 1;
    for(i = 0; i < gisecbuf.len; i++) {
      graph_load(&gisecbuf.b[i], gprefs, NULL);
      hash_table_print_stats(&db_graph.ht);
      graph_file_close(&gisecbuf.b[i]);
    }
    SWAP(db_graph.col_covgs, tmp_covgs);
    SWAP(db_graph.col_edges, isec_edges); db_graph.num_edge_cols = output_colours;
    // reset ginfo
    graph_info_init(&db_graph.ginfo[0]);
  }

  // Load graphs
  if(gfilebuf.len > 0)
  {
    GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);
    gprefs.must_exist_in_graph = (gisecbuf.len > 0);
    gprefs.must_exist_in_edges = isec_edges;

    for(i = 0; i < gfilebuf.len; i++) {
      graph_load(&gfilebuf.b[i], gprefs, NULL);
      hash_table_print_stats(&db_graph.ht);
      graph_file_close(&gfilebuf.b[i]);
    }
  }

  // Set sample names using seq_colours array
  for(i = 0; i < ncolours; i++) {
    strbuf_set(&db_graph.ginfo[samples[i].colour].sample_name, samples[i].name);
  }

  size_t start, end, num_load, colour, prev_colour = 0;

  // If we are using PCR duplicate removal,
  // it's best to load one colour at a time
  for(start = 0; start < ntasks; start = end, prev_colour = colour)
  {
    // Wipe read start bitfield
    colour = tasks[start].prefs.colour;
    if(remove_pcr_used)
    {
      if(colour != prev_colour)
        memset(db_graph.readstrt, 0, roundup_bits2bytes(db_graph.ht.capacity)*2);

      end = start+1;
      while(end < ntasks && end-start < MAX_IO_THREADS &&
            tasks[end].prefs.colour == colour) end++;
    }
    else {
      end = MIN2(start+MAX_IO_THREADS, ntasks);
    }

    num_load = end-start;
    build_graph(&db_graph, tasks+start, num_load, nthreads);
  }

  // Remove kmers with no coverage
  if(gisecbuf.len > 0) {
    db_graph_remove_no_covg_kmers(&db_graph, nthreads);
    db_graph_intersect_edges(&db_graph, nthreads, isec_edges);
  }

  // Print stats for hash table
  hash_table_print_stats(&db_graph.ht);

  // Print stats per input file
  for(i = 0; i < ntasks; i++) {
    build_graph_task_print_stats(&tasks[i]);
    build_graph_task_destroy(&tasks[i]);
  }

  status("Dumping graph...\n");
  graph_writer_save_mkhdr(out_path, &db_graph, sort_kmers, NULL, 0, output_colours);

  build_graph_task_buf_dealloc(&gtaskbuf);
  gfile_buf_dealloc(&gfilebuf);
  gfile_buf_dealloc(&gisecbuf);
  sample_name_buf_dealloc(&snamebuf);

  ctx_free(isec_edges);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
