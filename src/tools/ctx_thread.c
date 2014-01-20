#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "loading_stats.h"
#include "graph_info.h"
#include "graph_format.h"
#include "path_format.h"
#include "generate_paths.h"
#include "async_read_io.h"

#define DEFAULT_MIN_INS 0
#define DEFAULT_MAX_INS 500

static const char usage[] =
"usage: "CMD" thread [options] <out.ctp> <in.ctx>[:cols] [in2.ctx ...]\n"
"  Thread reads through the graph.  Save to file <out.ctp>.  <pop.ctx> can should\n"
"  have everyone in it and can be a pooled graph (with only 1 colour).  Samples\n"
"  are loaded from <in.ctx> files one at a time.\n"
"\n"
"  Options:\n"
"    -m <mem>                   How much memory to use\n"
"    -n <kmers>                 How many entries in the hash table\n"
"    -p <in.ctp>                Load existing path files first (multiple allowed)\n"
"    --col <colour>             Colour to thread through\n"
"    --seq <in.fa>              Thread reads from file (supports sam,bam,fq,*.gz)\n"
"    --seq2 <in.1.fq> <in.2.fq> Thread paired end reads\n"
"    --minIns <ins>             Minimum insert size for --seq2 [default:0]\n"
"    --maxIns <ins>             Maximum insert size for --seq2 [default:"QUOTE_MACRO(DEFAULT_MAX_INS)"]\n"
"    --twoway                   Use two-way gap filling\n"
"    --oneway                   Use one-way gap filling\n"
"\n"
"  When loading existing paths with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into.\n"
"\n"
"  Example:\n"
"    "CMD" thread -m 80G -p 0:Minnie.ctp -p 2:Mickey.ctp \\\n"
"                 --col 0 --seq2 sample3.1.fq sample3.2.fq \\\n"
"                 --col 1 --seq sample5.fa \\\n"
"                 samples.ctp samples.ctx\n"
"\n"
"  See `"CMD" pjoin` to combine .ctp files\n";

static void gen_paths_print_task(const GeneratePathsTask *t)
{
  int has_p2 = t->file2 != NULL;
  const char *p1 = t->file1->path, *p2 = has_p2 ? t->file2->path : "";
  char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

  if(t->fq_cutoff > 0) sprintf(fqCutoff, "%u", t->fq_cutoff);
  if(t->fq_offset > 0) sprintf(fqOffset, "%u", t->fq_offset);
  if(t->hp_cutoff > 0) sprintf(hpCutoff, "%u", t->hp_cutoff);

  status("[task] input: %s%s%s", p1, has_p2 ? ", " : "", p2);
  status("[task] insert min: %u max: %u", t->ins_gap_min, t->ins_gap_max);
  status("[task] FASTQ offset: %s, threshold: %s; cut homopolymers: %s",
         fqOffset, fqCutoff, hpCutoff);
  status("[task] %s-way gap traversal", t->one_way_gap_traverse ? "one" : "two");
  if(has_p2)
    status("[task] read pair: %s", t->read_pair_FR ? "FR" : "FF");
}

static void gen_path_task_create(const char *p1, const char *p2,
                                 size_t col, size_t min_ins, size_t max_ins,
                                 boolean one_way_bridge,
                                 uint32_t fq_offset, uint32_t fq_cutoff,
                                 uint32_t hp_cutoff, GeneratePathsTask *ptr)
{
  if(p1[0] == '-')
    print_usage(usage, "Path appears to be an option: %s", p1);
  if(p2 != NULL && p2[0] == '-')
    print_usage(usage, "Path appears to be an option: %s", p2);

  seq_file_t *f1, *f2 = NULL;

  if((f1 = seq_open(p1)) == NULL)
    die("Cannot read first --seq%s file: %s", p2 == NULL ? "" : "2", p1);
  if(p2 != NULL && (f2 = seq_open(p2)) == NULL)
    die("Cannot read second --seq2 file: %s", p2);

  GeneratePathsTask tsk = {.file1 = f1, .file2 = f2,
                           .ctxcol = 0, .ctpcol = col,
                           .ins_gap_min = min_ins, .ins_gap_max = max_ins,
                           .fq_offset = (uint8_t)fq_offset,
                           .fq_cutoff = (uint8_t)fq_cutoff,
                           .hp_cutoff = (uint8_t)hp_cutoff,
                           .read_pair_FR = true,
                           .one_way_gap_traverse = one_way_bridge};

  memcpy(ptr, &tsk, sizeof(GeneratePathsTask));
}

static int gen_path_task_cmp(const void *aa, const void *bb)
{
  const GeneratePathsTask *a = (const GeneratePathsTask*)aa;
  const GeneratePathsTask *b = (const GeneratePathsTask*)bb;
  if(a->ctpcol != b->ctpcol) return (int)a->ctpcol - (int)b->ctpcol;
  return (a > b ? 1 : (a < b ? -1 : 0));
}

static void gen_path_tasks_sort(GeneratePathsTask *tasks, size_t n)
{
  qsort(tasks, n, sizeof(GeneratePathsTask), gen_path_task_cmp);
}


static void get_binary_and_colour(const GraphFileReader *files, size_t num_graphs,
                                  size_t col, size_t *file_idx, size_t *col_idx)
{
  size_t i, n = 0;
  for(i = 0; i < num_graphs; i++) {
    if(n + graph_file_outncols(&files[i]) > col) {
      *col_idx = col - n; *file_idx = i; return;
    }
    n += graph_file_outncols(&files[i]);
  }
  die("Colour is greater than sum of graph colours [%zu > %zu]", col, n);
}

static int load_args(int argc, char **argv,
                     GeneratePathsTask *tasks, size_t *num_tasks_ptr)
{
  size_t num_tasks = 0;
  int argi;
  size_t min_ins = DEFAULT_MIN_INS, max_ins = DEFAULT_MAX_INS;
  uint32_t fq_offset = 0, fq_cutoff = 0, hp_cutoff = 0;
  boolean one_way_bridge = true;
  boolean col_set = false, col_used = false;
  size_t col;

  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcmp(argv[argi],"--printinsgaps") == 0)
    {
      // gen_paths_print_inserts is defined in add_read_paths.h
      gen_paths_print_inserts = true;
    }
    else if(strcmp(argv[argi],"--fq_threshold") == 0) {
      if(argi + 1 >= argc)
        print_usage(usage, "--fq_threshold <qual> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_offset) || fq_offset > 128)
        die("Invalid --fq_threshold argument: %s", argv[argi+1]);
      argi++;
    }
    else if(strcmp(argv[argi],"--fq_offset") == 0) {
      if(argi + 1 >= argc)
        print_usage(usage, "--fq_offset <offset> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &fq_offset) || fq_offset > 128)
        die("Invalid --fq_offset argument: %s", argv[argi+1]);
      argi++;
    }
    else if(strcmp(argv[argi],"--cut_hp") == 0) {
      if(argi + 1 >= argc)
        print_usage(usage, "--cut_hp <len> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &hp_cutoff))
        die("Invalid --cut_hp argument: %s", argv[argi+1]);
      if(hp_cutoff > UINT8_MAX)
        die("--cut_hp <hp> cannot be greater than %i", UINT8_MAX);
      argi++;
    }
    else if(strcmp(argv[argi],"--col") == 0)
    {
      if(argi+2 >= argc) print_usage(usage, "--col <colour> requires an argument");
      if(col_set && !col_used)
        print_usage(usage, "--seq or --seq2 must follow --col");
      if(!parse_entire_size(argv[argi+1], &col))
        print_usage(usage, "--col <colour> requires integers >= 0");
      col_set = true;
      col_used = false;
      argi++;
    }
    else if(strcasecmp(argv[argi],"--minIns") == 0)
    {
      if(argi+1 >= argc || !parse_entire_size(argv[++argi], &min_ins))
        print_usage(usage, "--minIns <bp> requires a positive integer argument");
      argi++;
    }
    else if(strcasecmp(argv[argi],"--maxIns") == 0)
    {
      if(argi+1 >= argc || !parse_entire_size(argv[++argi], &max_ins))
        print_usage(usage, "--maxIns <bp> requires a positive integer argument");
      if(max_ins < 20)
        warn("--maxGap < 20 seems very low!");
      argi++;
    }
    else if(strcasecmp(argv[argi],"--oneway") == 0) {
      one_way_bridge = true;
    }
    else if(strcasecmp(argv[argi],"--twoway") == 0) {
      one_way_bridge = false;
    }
    else if(strcmp(argv[argi],"--seq") == 0)
    {
      if(argi+1 == argc) print_usage(usage, "--seq <in.fa> missing args");
      if(!col_set) die("--seq <in.fa> before --col <colour>");

      gen_path_task_create(argv[argi+1], NULL, col,
                           min_ins, max_ins, one_way_bridge,
                           fq_offset, fq_cutoff, hp_cutoff,
                           &tasks[num_tasks]);

      num_tasks++;
      col_used = true;
      argi++;
    }
    else if(strcmp(argv[argi],"--seq2") == 0)
    {
      if(argi+2 >= argc) print_usage(usage, "--seq2 <in.1.fq> <in.2.fq> missing args");
      if(!col_set) die("--seq2 <in1.fa> <in2.fa> before --col <colour>");

      gen_path_task_create(argv[argi+1], argv[argi+2], col,
                           min_ins, max_ins, one_way_bridge,
                           fq_offset, fq_cutoff, hp_cutoff,
                           &tasks[num_tasks]);

      num_tasks++;
      col_used = true;
      argi += 2;
    }
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  if(num_tasks == 0)
    print_usage(usage, "need at least one: --col <c> --seq[2] <in> [in2]");
  if(!col_used)
    print_usage(usage, "--seq or --seq2 must follow last --col <col>");

  gen_path_tasks_sort(tasks, num_tasks);

  *num_tasks_ptr = num_tasks;
  return argi;
}

int ctx_thread(CmdArgs *args)
{
  cmd_accept_options(args, "atmnp", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  size_t max_tasks = (size_t)argc/2;
  GeneratePathsTask *tasks = malloc2(max_tasks * sizeof(GeneratePathsTask));
  size_t i, j, num_tasks, num_work_threads = args->max_work_threads;
  int argi; // arg index to continue from

  argi = load_args(argc, argv, tasks, &num_tasks);

  for(i = 0; i < num_tasks; i++)
    gen_paths_print_task(&tasks[i]);

  if(gen_paths_print_inserts && num_work_threads > 1) {
    die("--printinsgaps with >1 threads is a bad idea.");
  }

  if(argi + 1 >= argc) print_usage(usage, "Not enough arguments");

  char *out_ctp_path = argv[argi++];

  size_t num_graphs = (size_t)(argc - argi);
  char **graph_paths = argv + argi;

  //
  // Open graph graph_files
  //
  GraphFileReader graph_files[num_graphs];
  size_t ctx_max_kmers = 0, ctx_total_cols = 0;

  for(i = 0; i < num_graphs; i++)
  {
    status("File: %s", graph_paths[i]);
    graph_files[i] = INIT_GRAPH_READER;
    graph_file_open(&graph_files[i], graph_paths[i], true);

    if(i > 0 && graph_files[0].hdr.kmer_size != graph_files[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  graph_files[0].hdr.kmer_size, graph_files[i].hdr.kmer_size);
    }

    file_filter_update_intocol(&graph_files[i].fltr,
                               graph_files[i].fltr.intocol + ctx_total_cols);
    ctx_total_cols = graph_file_usedcols(&graph_files[i]);

    ctx_max_kmers = MAX2(ctx_max_kmers, graph_files[i].hdr.num_of_kmers);
  }

  //
  // Open path files
  //
  size_t num_pfiles = args->num_ctp_files;
  PathFileReader pfiles[num_pfiles];
  size_t path_max_mem = 0, path_max_usedcols = 0;

  for(i = 0; i < num_pfiles; i++) {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], args->ctp_files[i], true);
    path_max_mem = MAX2(path_max_mem, pfiles[i].hdr.num_path_bytes);
    path_max_usedcols = MAX2(path_max_usedcols, path_file_usedcols(&pfiles[i]));
  }

  // Check for path colours >= number of output path colours
  size_t max_gap_limit = 0;

  for(i = 0; i < num_tasks; i++) {
    max_gap_limit = MAX2(max_gap_limit, tasks[i].ins_gap_max);

    if(tasks[i].ctpcol >= ctx_total_cols) {
      print_usage(usage, "--col <C> >= number of path colours [%zu >= %zu]",
                  tasks[i].ctpcol, ctx_total_cols);
    }
  }

  size_t total_cols = MAX2(ctx_total_cols, path_max_usedcols);
  status("Creating paths file with %zu colours", total_cols);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem;
  char path_mem_str[100];

  bits_per_kmer = sizeof(Edges)*8 + sizeof(uint64_t)*8 +
                  2*num_work_threads + // Have traversed
                  1; // path store kmer lock

  // false -> don't use mem_to_use to decide how many kmers to store in hash
  // since we need some of that memory for storing paths
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, ctx_max_kmers,
                                        false, &graph_mem);

  // Path Memory
  size_t tmppathsize = paths_merge_needs_tmp(pfiles, num_pfiles) ? path_max_mem : 0;
  path_mem = args->mem_to_use - graph_mem;
  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s", path_mem_str);

  char req_path_mem_str[100];
  bytes_to_str(path_max_mem + tmppathsize, 1, req_path_mem_str);

  if(path_mem < path_max_mem + tmppathsize)
    die("Not enough memory for paths (requires %s)", req_path_mem_str);

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  //
  // Open output file
  //
  if(futil_file_exists(out_ctp_path))
    die("Output file already exists: %s", out_ctp_path);

  FILE *fout = fopen(out_ctp_path, "w");
  if(fout == NULL) die("Unable to open paths file to write: %s", out_ctp_path);
  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, graph_files[0].hdr.kmer_size, total_cols, 1, kmers_in_hash);
  kmers_in_hash = db_graph.ht.capacity;

  // path kmer locks
  db_graph.path_kmer_locks = calloc2(roundup_bits2bytes(kmers_in_hash), 1);

  // Edges
  db_graph.col_edges = calloc2(kmers_in_hash, sizeof(Edges));

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  PathStore *paths = &db_graph.pdata;
  uint8_t *path_store = malloc2(path_mem);
  path_store_init(paths, path_store, path_mem-tmppathsize, total_cols);
  uint8_t *path_tmpmem = path_store + path_mem - tmppathsize;

  // 1. Merge graph file headers into the graph
  size_t intocol, fromcol;
  for(i = 0; i < num_graphs; i++) {
    for(j = 0; j < graph_files[i].fltr.ncols; j++) {
      intocol = graph_file_intocol(&graph_files[i], j);
      fromcol = graph_file_fromcol(&graph_files[i], j);
      graph_info_merge(&db_graph.ginfo[intocol], &graph_files[i].hdr.ginfo[fromcol]);
    }
  }

  // Load paths
  if(num_pfiles > 0) {
    // Paths loaded into empty colours will update the sample names
    paths_format_merge(pfiles, num_pfiles, true, path_tmpmem, tmppathsize, &db_graph);
    path_store_resize(paths, path_mem);
  }

  // Set up paths header. This is for the output file we are creating
  PathFileHeader pheader = INIT_PATH_FILE_HDR_MACRO;
  paths_header_alloc(&pheader, total_cols);

  pheader.num_of_cols = (uint32_t)total_cols;
  pheader.kmer_size = graph_files[0].hdr.kmer_size;

  // Set new path header sample names from graph header
  for(i = 0; i < total_cols; i++)
    strbuf_set(&pheader.sample_names[i], db_graph.ginfo[i].sample_name.buff);

  // 2. reduce number of graph colours
  db_graph_realloc(&db_graph, 1, 1);

  // Setup for loading graphs graph
  SeqLoadingStats *gstats = seq_loading_stats_create(0);
  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = false};

  //
  // Start up the threads, do the work
  //
  GenPathWorker *workers;
  workers = gen_paths_workers_alloc(num_work_threads, &db_graph);

  // ... Send jobs ...
  size_t start, end;
  for(start = 0; start < num_tasks; start = end)
  {
    // wipe colour 0
    graph_info_init(&db_graph.ginfo[0]);
    memset(db_graph.col_edges, 0, kmers_in_hash * sizeof(Edges));

    // Load graph
    size_t ctxindex, ctxcol;
    size_t ctpcol = tasks[start].ctpcol;

    get_binary_and_colour(graph_files, num_graphs, ctpcol, &ctxindex, &ctxcol);
    graph_load_colour(&graph_files[ctxindex], gprefs, gstats, ctxcol, 0);

    // Get list of input files to read
    end = start+1;
    while(end < num_tasks && end-start < args->max_io_threads &&
          tasks[end].ctpcol == ctpcol) end++;

    generate_paths(tasks+start, end-start, workers, num_work_threads);
  }

  // Output statistics
  const uint64_t *ins_gaps, *err_gaps;
  size_t ins_gaps_len, err_gaps_len;

  ins_gaps = gen_paths_get_ins_gap(workers, &ins_gaps_len);
  err_gaps = gen_paths_get_err_gap(workers, &err_gaps_len);

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  gen_paths_get_stats(workers, num_work_threads, stats);

  // Print mp gap size / insert stats to a file
  gen_paths_dump_gap_sizes("gap_sizes.%u.csv", err_gaps, err_gaps_len,
                           db_graph.kmer_size, false,
                           stats->num_se_reads + stats->num_pe_reads);

  if(stats->num_pe_reads > 0) {
    gen_paths_dump_gap_sizes("mp_sizes.%u.csv", ins_gaps, ins_gaps_len,
                             db_graph.kmer_size, true, stats->num_pe_reads);
  }
  else
    status("No PE reads parsed");

  // ins_gap, err_gap no longer allocated after this line
  gen_paths_workers_dealloc(workers, num_work_threads);

  // Done

  char se_num_str[100], pe_num_str[100], sepe_num_str[100];
  ulong_to_str(stats->num_se_reads, se_num_str);
  ulong_to_str(stats->num_pe_reads / 2, pe_num_str);
  ulong_to_str(stats->num_se_reads + stats->num_pe_reads, sepe_num_str);
  status("[stats] single reads: %s; read pairs: %s; total: %s",
         se_num_str, pe_num_str, sepe_num_str);

  size_t num_path_bytes = (size_t)(paths->next - paths->store);
  char kmers_str[100], paths_str[100], mem_str[100];
  ulong_to_str(paths->num_kmers_with_paths, kmers_str);
  ulong_to_str(paths->num_of_paths, paths_str);
  bytes_to_str(num_path_bytes, 1, mem_str);

  status("Saving paths: %s paths, %s path-bytes, %s kmers",
         paths_str, mem_str, kmers_str);

  // Update header and write
  paths_header_update(&pheader, paths);
  paths_format_write_header(&pheader, fout);
  paths_format_write_optimised_paths(&db_graph, fout);

  fclose(fout);

  free(tasks);

  free(db_graph.col_edges);
  free((void *)db_graph.kmer_paths);
  free((void *)db_graph.path_kmer_locks);
  free(path_store);

  paths_header_dealloc(&pheader);

  seq_loading_stats_free(stats);
  seq_loading_stats_free(gstats);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_graphs; i++) graph_file_dealloc(&graph_files[i]);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  status("Paths written to: %s", out_ctp_path);

  return EXIT_SUCCESS;
}
