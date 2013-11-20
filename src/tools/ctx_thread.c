#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "add_read_paths.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_walker.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" thread [options] <out.ctp> <in.ctx>[:cols] [in2.ctx ...]\n"
"  Thread reads through the graph.  Save to file <out.ctp>.  <pop.ctx> can should\n"
"  have everyone in it and can be a pooled graph (with only 1 colour).  Samples\n"
"  are loaded from <in.ctx> files one at a time.\n"
"\n"
"  Options:\n"
"    -m <mem>                   How much memory to use\n"
"    -n <kmers>                 How many entries in the hash table\n"
"    --col <colour>             Colour to thread through\n"
"    --seq <in.fa>              Thread reads from file (supports sam,bam,fq,*.gz)\n"
"    --seq2 <in.1.fq> <in.2.fq> Thread paired end reads\n"
"\n"
// "  Insert size filtering for follow --seq2 files:\n"
// "    --minIns <ins>             Minimum insert size [default:0]\n"
// "    --maxIns <ins>             Maximum insert size [default:500]\n"
// "\n"
"  Example: If we want paths for 2 samples only we can do:\n"
"    "CMD" thread -m 80G \\\n"
"                 --col 0 --seq2 sample3.1.fq sample3.2.fq \\\n"
"                 --col 1 --seq sample5.fa \\\n"
"                 sample3and5.ctp population.c6.ctx:2 population.c6.ctx:4\n"
"\n"
"  Or you could pool the samples first:\n"
"    "CMD" join --flatten pool.samples5and3.ctx population.c6.ctx:3,5\n"
"    "CMD" thread -m 80G --col 1 --seq2 sample3.1.fq sample3.2.fq \\\n"
"                 2 sample3.3and5.ctp pool.samples5and3.ctx\n"
"    "CMD" thread -m 80G --col 0 --seq sample5.fa \\\n"
"                 2 sample5.3and5.ctp pool.samples5and3.ctx population.c6.ctx\n"
"    "CMD" pmerge sample3and5.ctp sample3.3and5.ctp sample5.3and5.ctp\n";

#define NUM_PASSES 1

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
  die("Colour is greater than sum of binary colours [%zu > %zu]", col, n);
}

int ctx_thread(CmdArgs *args)
{
  cmd_accept_options(args, "tmn", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  size_t graph_col, num_seq_cols = 0, graph_colours[argc];

  seq_file_t *seqfiles[argc];
  size_t num_sf = 0, sf = 0;

  int argi;
  boolean used_last_path = false;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcasecmp(argv[argi],"--col") == 0)
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
    else if(strcasecmp(argv[argi],"--seq") == 0)
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
    else if(strcasecmp(argv[argi],"--seq2") == 0)
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
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  int argend = argi;

  if(argend + 1 >= argc) print_usage(usage, "Not enough arguments");

  char *out_ctp_path = argv[argi++];

  if(file_exists(out_ctp_path))
    die("Output file already exists: %s", out_ctp_path);

  size_t num_files = argc - argi;
  char **graph_paths = argv + argi;

  //
  // Probe graph files
  //
  GraphFileReader files[num_files];
  size_t i, j, ctx_max_kmers = 0, total_cols = 0;

  // Set up paths header
  PathFileHeader pheader = {.version = CTX_PATH_FILEFORMAT,
                            .num_of_cols = 0,
                            .capacity = 0};

  // Validate input graphs
  for(i = 0; i < num_files; i++)
  {
    status("File: %s", graph_paths[i]);
    files[i] = INIT_GRAPH_READER;
    int ret = graph_file_open(&files[i], graph_paths[i], false);

    if(ret == 0)
      print_usage(usage, "Cannot read input graph file: %s", graph_paths[i]);
    else if(ret < 0)
      print_usage(usage, "Input graph file isn't valid: %s", graph_paths[i]);

    if(i > 0 && files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, files[i].hdr.kmer_size);
    }

    size_t ncols = files[i].intocol + graph_file_outncols(&files[i]);
    files[i].intocol += total_cols;
    total_cols += ncols;

    ctx_max_kmers = MAX2(ctx_max_kmers, files[i].hdr.num_of_kmers);
  }

  status("Total %zu cols", total_cols);

  pheader.kmer_size = files[0].hdr.kmer_size;
  paths_header_alloc(&pheader, total_cols);
  pheader.num_of_cols = total_cols;

  // Set path header sample names
  for(i = 0; i < num_files; i++) {
    for(j = 0; j < files[i].ncols; j++) {
      strbuf_set(&pheader.sample_names[files[i].intocol+j],
                 files[i].hdr.ginfo[files[i].cols[j]].sample_name.buff);
    }
  }

  // Check for invalid or duplicate ctp colours
  // or path colours >= number of output path colours
  qsort(graph_colours, num_seq_cols, sizeof(size_t), cmp_size);
  for(i = 0; i < num_seq_cols; i++) {
    if(i+1 < num_seq_cols && graph_colours[i] == graph_colours[i+1])
      print_usage(usage, "Duplicate --col <ctpcol> given: %zu", graph_colours[i]);
    if(graph_colours[i] >= pheader.num_of_cols) {
      print_usage(usage, "output path colour >= number of path colours [%zu >= %u]",
                  graph_colours[i], pheader.num_of_cols);
    }
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem;
  char graph_mem_str[100], path_mem_str[100];

  bits_per_kmer = sizeof(Edges)*8 + sizeof(uint64_t)*8 +
                  2*args->num_threads; // Have traversed

  // false -> don't use mem_to_use to decide how many kmers to store in hash
  // since we need some of that memory for storing paths
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, ctx_max_kmers, false);

  graph_mem = hash_table_mem(kmers_in_hash,false,NULL) +
              (kmers_in_hash*bits_per_kmer)/8;
  bytes_to_str(graph_mem, 1, graph_mem_str);

  if(graph_mem >= args->mem_to_use)
    die("Not enough memory for graph (requires %s)", graph_mem_str);

  // Path Memory
  path_mem = args->mem_to_use - graph_mem;
  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  // Open output file
  FILE *fout = fopen(out_ctp_path, "w");

  if(fout == NULL)
    die("Unable to open binary paths file to write: %s\n", out_ctp_path);

  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, pheader.kmer_size, 1, 1, kmers_in_hash);
  kmers_in_hash = db_graph.ht.capacity;

  // Edges
  db_graph.col_edges = calloc2(kmers_in_hash, sizeof(uint8_t));

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  path_store_init(&db_graph.pdata, path_store, path_mem, pheader.num_of_cols);

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

  uint32_t gap_limit = 500;
  PathsWorkerPool *pool;

  pool = paths_worker_pool_new(args->num_threads, &db_graph, gap_limit);

  // Parse input sequence
  status("Threading reads through the graph...\n");

  size_t rep, ctxindex, ctxcol;

  for(rep = 0; rep < NUM_PASSES; rep++)
  {
    for(argi = 0; argi < argend; argi++)
    {
      if(strcmp(argv[argi], "--col") == 0)
      {
        // wipe colour 0
        // db_graph_wipe_colour(&db_graph, 0);
        graph_info_init(&db_graph.ginfo[0]);
        memset(db_graph.col_edges, 0, kmers_in_hash * sizeof(uint8_t));

        parse_entire_size(argv[argi+1], &graph_col);
        add_paths_set_colours(pool, 0, graph_col); // no need to pass second col

        // Pick correct binary and colour
        get_binary_and_colour(files, num_files, graph_col, &ctxindex, &ctxcol);
        status("Load file with index %zu and colour index %zu [%zu]",
               ctxindex, ctxcol, graph_col);
        graph_load_colour(&files[ctxindex], &prefs, stats, ctxcol, 0);

        argi++;
      }
      else if(strcmp(argv[argi], "--seq") == 0) {
        add_read_paths_to_graph(pool, seqfiles[sf], NULL, gap_limit,
                                &prefs, stats);
        argi += 1;
        sf++;
      }
      else if(strcmp(argv[argi], "--seq2") == 0) {
        add_read_paths_to_graph(pool, seqfiles[sf], seqfiles[sf+1], gap_limit,
                                &prefs, stats);
        argi += 2;
        sf += 2;
      }
      else die("Unknown arg: %s", argv[argi]);
    }
  }

  paths_worker_pool_dealloc(pool);

  PathStore *paths = &db_graph.pdata;
  size_t num_path_bytes = paths->next - paths->store;
  char kmers_str[100], paths_str[100], mem_str[100];
  ulong_to_str(paths->num_kmers_with_paths, kmers_str);
  ulong_to_str(paths->num_of_paths, paths_str);
  bytes_to_str(num_path_bytes, 1, mem_str);

  status("Saving paths: %s paths, %s path-bytes, %s kmers\n",
          paths_str, mem_str, kmers_str);

  paths_format_write_optimised_paths(&db_graph, fout);

  // Update header and overwrite
  paths_header_update(&pheader, &db_graph.pdata);
  fseek(fout, 0, SEEK_SET);
  paths_format_write_header_core(&pheader, fout);

  fclose(fout);

  free(db_graph.col_edges);
  free((void *)db_graph.kmer_paths);
  free(path_store);

  paths_header_dealloc(&pheader);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);

  status("Paths written to: %s\n", out_ctp_path);

  return EXIT_SUCCESS;
}
