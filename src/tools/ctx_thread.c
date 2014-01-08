#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "add_read_paths.h"
#include "add_path_workers.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_walker.h"
#include "seq_reader.h"

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
"\n"
"  When loading existing paths with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into.\n"
"\n"
"  Insert size filtering for following --seq2 files:\n"
"    --minIns <ins>             Minimum insert size [default:0]\n"
"    --maxIns <ins>             Maximum insert size [default:"QUOTE_MACRO(DEFAULT_MAX_INS)"]\n"
"\n"
"  Example:\n"
"    "CMD" thread -m 80G -p 0:Minnie.ctp -p 2:Mickey.ctp \\\n"
"                 --col 0 --seq2 sample3.1.fq sample3.2.fq \\\n"
"                 --col 1 --seq sample5.fa \\\n"
"                 samples.ctp samples.ctx\n"
"\n"
"  See `"CMD" pjoin` to combine .ctp files\n";

static void get_binary_and_colour(const GraphFileReader *files, size_t num_files,
                                  size_t col, size_t *file_idx, size_t *col_idx)
{
  size_t i, n = 0;
  for(i = 0; i < num_files; i++) {
    if(n + graph_file_outncols(&files[i]) > col) {
      *col_idx = col - n; *file_idx = i; return;
    }
    n += graph_file_outncols(&files[i]);
  }
  die("Colour is greater than sum of graph colours [%zu > %zu]", col, n);
}

int ctx_thread(CmdArgs *args)
{
  cmd_accept_options(args, "tmnp", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  size_t graph_col, num_seq_cols = 0, graph_colours[argc];

  seq_file_t *seqfiles[argc];
  size_t num_sf = 0, sf = 0;
  size_t gap_limit = DEFAULT_MAX_INS, gap_tmp;

  int argi;
  boolean used_last_path = false;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcmp(argv[argi],"--printinsgaps") == 0)
    {
      // print_traversed_inserts is defined in add_read_paths.h
      print_traversed_inserts = true;
    }
    else if(strcmp(argv[argi],"--col") == 0)
    {
      if(argi+2 >= argc)
        print_usage(usage, "--col <ctxcol> <ctpcol> requires an argument");
      if(num_seq_cols > 0 && !used_last_path)
        print_usage(usage, "--seq or --seq2 must follow --col");
      if(!parse_entire_size(argv[argi+1], &graph_col))
        print_usage(usage, "--col <ctxcol> requires integers >= 0");
      graph_colours[num_seq_cols++] = graph_col;
      used_last_path = false;
      argi++;
    }
    else if(strcmp(argv[argi],"--seq") == 0)
    {
      if(argi+1 == argc)
        print_usage(usage, "--seq <in.fa> requires an argument");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read --seq file: %s", argv[argi+1]);
      if(num_seq_cols == 0)
        die("--seq <in.fa> before --col <ctxcol> <ctpcol>");
      used_last_path = true;
      argi++;
    }
    else if(strcmp(argv[argi],"--seq2") == 0)
    {
      if(argi+2 >= argc)
        print_usage(usage, "--seq2 <in.1.fq> <in.2.fq> requires two arguments");
      if(num_seq_cols == 0)
        die("--seq2 <in1.fa> <in2.fa> before --col <ctxcol> <ctpcol>");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read first --seq2 file: %s", argv[argi+1]);
      if((seqfiles[num_sf++] = seq_open(argv[argi+2])) == NULL)
        die("Cannot read second --seq2 file: %s", argv[argi+2]);
      used_last_path = true;
      argi += 2;
    }
    else if(strcasecmp(argv[argi],"--minIns") == 0)
    {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &gap_tmp))
        print_usage(usage, "--minIns <bp> requires a positive integer argument");
      argi++;
    }
    else if(strcasecmp(argv[argi],"--maxIns") == 0)
    {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &gap_tmp))
        print_usage(usage, "--maxIns <bp> requires a positive integer argument");
      gap_limit = MAX2(gap_limit, gap_tmp);
      argi++;
    }
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  if(print_traversed_inserts && args->num_threads > 1) {
    die("--printinsgaps with >1 threads is a bad idea.");
  }

  if(num_seq_cols == 0)
    print_usage(usage, "need at least one: --col <c> --seq[2] <in> [in2]");
  if(num_seq_cols > 0 && !used_last_path)
    print_usage(usage, "--seq or --seq2 must follow --col");

  int argend = argi;

  if(argend + 1 >= argc) print_usage(usage, "Not enough arguments");

  char *out_ctp_path = argv[argi++];

  if(futil_file_exists(out_ctp_path))
    die("Output file already exists: %s", out_ctp_path);

  size_t num_files = (size_t)(argc - argi);
  char **graph_paths = argv + argi;

  //
  // Open graph files
  //
  GraphFileReader files[num_files];
  size_t i, j, ctx_max_kmers = 0, total_cols = 0;

  for(i = 0; i < num_files; i++)
  {
    status("File: %s", graph_paths[i]);
    files[i] = INIT_GRAPH_READER;
    graph_file_open(&files[i], graph_paths[i], true);

    if(i > 0 && files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, files[i].hdr.kmer_size);
    }

    file_filter_update_intocol(&files[i].fltr, files[i].fltr.intocol + total_cols);
    total_cols = graph_file_usedcols(&files[i]);

    ctx_max_kmers = MAX2(ctx_max_kmers, files[i].hdr.num_of_kmers);
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

  total_cols = MAX2(total_cols, path_max_usedcols);
  status("Creating paths file with %zu colours", total_cols);

  // Check for invalid or duplicate ctp colours
  // or path colours >= number of output path colours
  qsort(graph_colours, num_seq_cols, sizeof(size_t), cmp_size);
  for(i = 0; i < num_seq_cols; i++) {
    if(i+1 < num_seq_cols && graph_colours[i] == graph_colours[i+1])
      print_usage(usage, "Duplicate --col <ctpcol> given: %zu", graph_colours[i]);
    if(graph_colours[i] >= total_cols) {
      print_usage(usage, "output path colour >= number of path colours [%zu >= %zu]",
                  graph_colours[i], total_cols);
    }
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem;
  char path_mem_str[100];

  bits_per_kmer = sizeof(Edges)*8 + sizeof(uint64_t)*8 +
                  2*args->num_threads; // Have traversed

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

  // Open output file
  FILE *fout = fopen(out_ctp_path, "w");

  if(fout == NULL)
    die("Unable to open paths file to write: %s\n", out_ctp_path);

  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, total_cols, 1, kmers_in_hash);
  kmers_in_hash = db_graph.ht.capacity;

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
  for(i = 0; i < num_files; i++) {
    for(j = 0; j < files[i].fltr.ncols; j++) {
      intocol = graph_file_intocol(&files[i], j);
      fromcol = graph_file_fromcol(&files[i], j);
      graph_info_merge(&db_graph.ginfo[intocol], &files[i].hdr.ginfo[fromcol]);
    }
  }

  // Load paths
  if(num_pfiles > 0) {
    // Paths loaded into empty colours will update the sample names
    paths_format_merge(pfiles, num_pfiles, true, path_tmpmem, tmppathsize, &db_graph);
    path_store_resize(paths, path_mem);
  }

  // Set up paths header. This is for the output file we are creating
  PathFileHeader pheader = {.version = CTX_PATH_FILEFORMAT,
                            .kmer_size = files[0].hdr.kmer_size,
                            .num_of_cols = 0,
                            .capacity = 0};

  paths_header_alloc(&pheader, total_cols);
  pheader.num_of_cols = (uint32_t)total_cols;

  // Set new path header sample names from graph header
  for(i = 0; i < total_cols; i++)
    strbuf_set(&pheader.sample_names[i], db_graph.ginfo[i].sample_name.buff);

  // 2. reduce number of graph colours
  db_graph_realloc(&db_graph, 1, 1);

  // Setup for loading graphs graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           // binaries
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = false,
                           // Sequence
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false};

  paths_format_write_header(&pheader, fout);

  PathsWorkerPool *pool;

  // set up mutexes for reading paths
  add_read_paths_init();
  pool = path_workers_pool_new(args->num_threads, &db_graph, gap_limit+1);

  // Parse input sequence
  size_t ctxindex, ctxcol, min_ins = 0, max_ins = DEFAULT_MAX_INS;

  for(argi = 0; argi < argend; argi++)
  {
    if(strcmp(argv[argi],"--printinsgaps")) {}
    else if(strcmp(argv[argi], "--col") == 0)
    {
      if(argi > 0) path_workers_wait_til_finished(pool);

      // wipe colour 0
      graph_info_init(&db_graph.ginfo[0]);
      memset(db_graph.col_edges, 0, kmers_in_hash * sizeof(Edges));

      parse_entire_size(argv[argi+1], &graph_col);

      // Pick correct graph file and colour
      get_binary_and_colour(files, num_files, graph_col, &ctxindex, &ctxcol);
      graph_load_colour(&files[ctxindex], &prefs, stats, ctxcol, 0);

      argi++;
    }
    else if(strcasecmp(argv[argi],"--minIns") == 0)
    {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &min_ins))
        die("Bad --minIns");
      argi++;
    }
    else if(strcasecmp(argv[argi],"--maxIns") == 0)
    {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &max_ins))
        die("Bad --maxIns");
      argi++;
    }
    else if(strcmp(argv[argi], "--seq") == 0) {
      path_workers_add_paths_to_graph(pool, seqfiles[sf], NULL,
                                      0, graph_col, min_ins, max_ins,
                                      &prefs, stats);
      argi += 1;
      sf++;
    }
    else if(strcmp(argv[argi], "--seq2") == 0) {
      path_workers_add_paths_to_graph(pool, seqfiles[sf], seqfiles[sf+1],
                                      0, graph_col, min_ins, max_ins,
                                      &prefs, stats);
      argi += 2;
      sf += 2;
    }
    else die("Unknown arg: %s", argv[argi]);
  }

  path_workers_pool_dealloc(pool, stats->num_se_reads, stats->num_pe_reads / 2);

  char se_num_str[100], pe_num_str[100], sepe_num_str[100];
  ulong_to_str(stats->num_se_reads, se_num_str);
  ulong_to_str(stats->num_pe_reads / 2, pe_num_str);
  ulong_to_str(stats->num_se_reads + stats->num_pe_reads, sepe_num_str);
  status("Threaded: single reads: %s; read pairs: %s; total: %s",
         se_num_str, pe_num_str, sepe_num_str);

  size_t num_path_bytes = (size_t)(paths->next - paths->store);
  char kmers_str[100], paths_str[100], mem_str[100];
  ulong_to_str(paths->num_kmers_with_paths, kmers_str);
  ulong_to_str(paths->num_of_paths, paths_str);
  bytes_to_str(num_path_bytes, 1, mem_str);

  status("Saving paths: %s paths, %s path-bytes, %s kmers",
         paths_str, mem_str, kmers_str);

  paths_format_write_optimised_paths(&db_graph, fout);

  // Update header and overwrite
  paths_header_update(&pheader, paths);
  fseek(fout, 0, SEEK_SET);
  paths_format_write_header_core(&pheader, fout);

  fclose(fout);

  add_read_paths_cleanup();

  free(db_graph.col_edges);
  free((void *)db_graph.kmer_paths);
  free(path_store);

  paths_header_dealloc(&pheader);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  status("Paths written to: %s", out_ctp_path);

  return EXIT_SUCCESS;
}
