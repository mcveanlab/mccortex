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
"    -h <kmers>                 How many entries in the hash table\n"
"    --col <ctxcol> <ctpcol>    Load ctx path into ctp col\n"
"    --seq <in.fa>              Thread reads from file (supports sam,bam,fq,*.gz)\n"
"    --seq2 <in.1.fq> <in.2.fq> Thread paired end reads\n"
"\n"
// "  Insert size filtering for follow --seq2 files:\n"
// "    --minIns <ins>             Minimum insert size [default:0]\n"
// "    --maxIns <ins>             Maximum insert size [default:500]\n"
// "\n"
"  Example: If we want paths for 2 samples only we can do:\n"
"    "CMD" thread -m 80G \\\n"
"                 --col 3 0 --seq2 sample3.1.fq sample3.2.fq \\\n"
"                 --col 5 1 --seq sample5.fa \\\n"
"                 2 sample3and5.ctp population.c6.ctx population.c6.ctx\n"
"\n"
"  Or you could pool the samples first:\n"
"    "CMD" join --flatten pool.samples5and3.ctx population.c6.ctx:3,5\n"
"    "CMD" thread -m 80G --col 3 0 --seq2 sample3.1.fq sample3.2.fq \\\n"
"                 2 sample3.3and5.ctp pool.samples5and3.ctx population.c6.ctx\n"
"    "CMD" thread -m 80G --col 5 1 --seq sample5.fa \\\n"
"                 2 sample5.3and5.ctp pool.samples5and3.ctx population.c6.ctx\n"
"    "CMD" pmerge sample3and5.ctp sample3.3and5.ctp sample5.3and5.ctp\n";

#define NUM_PASSES 1

static void get_binary_and_colour(uint32_t colour, uint32_t num_binaries,
                                  uint32_t ctx_cols[num_binaries],
                                  uint32_t *ctxindex, uint32_t *ctxcol)
{
  uint32_t i, colsum = 0;
  for(i = 0; i < num_binaries; i++)
  {
    colsum += ctx_cols[i];
    if(colour < colsum) {
      *ctxindex = i;
      *ctxcol = colour - (colsum - ctx_cols[i]);
      return;
    }
  }
  die("Colour is greater than sum of binary colours [%u > %u]", colour, colsum);
}

int ctx_thread(CmdArgs *args)
{
  cmd_accept_options(args, "tm");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  uint32_t graph_col, path_col, num_of_path_cols = 0;
  uint32_t path_colours[argc], graph_colours[argc];

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
      if(num_of_path_cols > 0 && !used_last_path)
        print_usage(usage, "--seq or --seq2 must follow --col");
      if(!parse_entire_uint(argv[argi+1], &graph_col))
        print_usage(usage, "--col <ctxcol> <ctpcol> requires integers >= 0");
      if(!parse_entire_uint(argv[argi+2], &path_col))
        print_usage(usage, "--col <ctxcol> <ctpcol> requires integers >= 0");
      graph_colours[num_of_path_cols] = graph_col;
      path_colours[num_of_path_cols++] = path_col;
      used_last_path = false;
      argi += 2;
    }
    else if(strcasecmp(argv[argi],"--seq") == 0)
    {
      if(argi+1 == argc)
        print_usage(usage, "--seq <in.fa> requires an argument");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read --seq file: %s", argv[argi+1]);
      if(num_of_path_cols == 0)
        die("--seq <in.fa> before --col <ctxcol> <ctpcol>");
      used_last_path = true;
      argi++;
    }
    else if(strcasecmp(argv[argi],"--seq2") == 0)
    {
      if(argi+2 >= argc)
        print_usage(usage, "--seq2 <in.1.fq> <in.2.fq> requires two arguments");
      if(num_of_path_cols == 0)
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

  if(argend + 2 > argc) print_usage(usage, "Not enough arguments");

  char *out_ctp_path = argv[argi++];

  if(file_exists(out_ctp_path))
    die("Output file already exists: %s", out_ctp_path);

  size_t num_binaries = argc - argi;
  char **binary_paths = argv + argi;

  //
  // Probe graph files
  //
  boolean is_binary = false;
  uint32_t ctx_num_of_cols[num_binaries], ctx_max_cols[num_binaries];
  uint32_t max_num_cols = 0;

  GraphFileHeader gheader = INIT_GRAPH_FILE_HDR;

  uint64_t ctx_max_kmers = 0;
  size_t i;

  // Set up paths header
  PathFileHeader pheader = {.version = CTX_PATH_FILEFORMAT,
                            .num_of_cols = 0,
                            .capacity = 0};

  // Validate input graphs
  for(i = 0; i < num_binaries; i++)
  {
    if(!graph_file_probe(binary_paths[i], &is_binary, &gheader))
      print_usage(usage, "Cannot read input graph file: %s", binary_paths[i]);
    else if(!is_binary)
      print_usage(usage, "Input graph file isn't valid: %s", binary_paths[i]);

    if(i > 0 && pheader.kmer_size != gheader.kmer_size) {
      die("Graph kmer-sizes do not match [%u vs %u; %s]\n",
          pheader.kmer_size, gheader.kmer_size, binary_paths[i]);
    }
    pheader.kmer_size = gheader.kmer_size;

    ctx_num_of_cols[i] = gheader.num_of_cols;
    ctx_max_cols[i] = gheader.max_col;
    max_num_cols = MAX2(gheader.num_of_cols, max_num_cols);
    ctx_max_kmers = MAX2(ctx_max_kmers, gheader.num_of_kmers);

    // alloc also reallocs
    paths_header_alloc(&pheader, pheader.num_of_cols + gheader.num_of_cols);

    // Clumsy way to update sample names
    uint32_t col;
    for(col = 0; col < gheader.num_of_cols; col++)
    {
      // gheader.ginfo[col] should be gheader.ginfo[filter.cols[col]]
      strbuf_set(&pheader.sample_names[pheader.num_of_cols+col],
                 gheader.ginfo[col].sample_name.buff);
    }

    pheader.num_of_cols += gheader.num_of_cols;
  }

  // Get colour indices
  uint32_t load_colours[num_binaries][max_num_cols];
  for(i = 0; i < num_binaries; i++)
    graph_file_parse_colours(binary_paths[i], load_colours[i], ctx_max_cols[i]);


  // Check for invalid or duplicate ctp colours
  // or path colours >= number of output path colours
  qsort(path_colours, num_of_path_cols, sizeof(uint32_t), cmp_uint32);
  for(i = 0; i < num_of_path_cols; i++) {
    if(i+1 < num_of_path_cols && path_colours[i] == path_colours[i+1])
      print_usage(usage, "Duplicate --col <ctpcol> given: %u", path_colours[i]);
    if(path_colours[i] >= pheader.num_of_cols) {
      print_usage(usage, "output path colour >= number of path colours [%u >= %u]",
                  path_colours[i], pheader.num_of_cols);
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

  // Edges
  db_graph.col_edges = calloc2(kmers_in_hash, sizeof(uint8_t));

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  path_store_init(&db_graph.pdata, path_store, path_mem, pheader.num_of_cols);

  // In colour (only need for one colour)
  // size_t words64_per_col = round_bits_to_words64(kmers_in_hash);
  // db_graph.node_in_cols = calloc2(cols_in_graph * words64_per_col, sizeof(uint64_t));

  // Setup for loading graphs graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .db_graph = &db_graph,
                           // binaries
                           .merge_colours = true,
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

  size_t rep;
  uint32_t ctxindex, ctxcol;

  for(rep = 0; rep < NUM_PASSES; rep++)
  {
    for(argi = 0; argi < argend; argi++)
    {
      if(strcmp(argv[argi], "--col") == 0)
      {
        // wipe colour 0
        db_graph_wipe_colour(&db_graph, 0);
        graph_info_init(&db_graph.ginfo[0]);

        parse_entire_uint(argv[argi+1], &graph_col);
        parse_entire_uint(argv[argi+2], &path_col);
        add_paths_set_colours(pool, 0, path_col); // no need to pass second col

        // Pick correct binary and colour
        get_binary_and_colour(graph_col, num_binaries, ctx_num_of_cols,
                              &ctxindex, &ctxcol);

        graph_load_colour(binary_paths[ctxindex], &prefs, stats,
                          load_colours[ctxindex][ctxcol]);

        strbuf_set(&pheader.sample_names[path_col],
                   db_graph.ginfo[0].sample_name.buff);

        argi += 2;
      }
      else if(strcmp(argv[argi], "--seq") == 0) {
        add_read_paths_to_graph(pool, seqfiles[sf++], NULL, gap_limit,
                                &prefs, stats);
        argi += 1;
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
  // free(db_graph.node_in_cols);
  free((void *)db_graph.kmer_paths);
  free(path_store);

  graph_header_dealloc(&gheader);
  paths_header_dealloc(&pheader);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  status("Paths written to: %s\n", out_ctp_path);

  return EXIT_SUCCESS;
}
