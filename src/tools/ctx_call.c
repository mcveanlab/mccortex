#include "global.h"

#include "tools.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_paths.h"
#include "graph_walker.h"
#include "bubble_caller.h"
#include "db_node.h"

// DEV: "usage: "CMD" call [options] <out.bubbles.gz> <in.ctx> [in2.ctx ...]\n"
const char call_usage[] =
"usage: "CMD" call [options] <in.ctx> <out.bubbles.gz>\n"
"  Find bubbles (potential variants) in graph file in.ctx, save to out.bubbles.gz\n"
"\n"
"  Options:\n"
"    -m <memory> | -t <threads> | -p <paths.ctp>\n"
"    -p <in.ctp>        Load path file (can specify multiple times)\n"
"    --ref <col>        Reference genome in given colour (can specify multiple times)\n"
"    --maxallele <len>  Max bubble branch length in kmers [default: 300]\n"
"    --maxflank <len>   Max flank length in kmers [default: 1000]\n"
"\n"
"  When loading path files with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into.\n";

int ctx_call(CmdArgs *args)
{
  int argi, argc = args->argc;
  char **argv = args->argv;
  if(argc < 2 || argc & 1) cmd_print_usage(NULL);

  size_t num_of_threads = args->max_work_threads;
  size_t i, ref_cols[argc], num_ref = 0;
  size_t max_allele_len = 300, max_flank_len = 1000;

  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcmp(argv[argi],"--ref") == 0) {
      if(argi + 1 == argc || !parse_entire_size(argv[argi+1], &ref_cols[num_ref]))
        cmd_print_usage("--ref <col> requires an int arg");
      num_ref++; argi++;
    }
    else if(strcmp(argv[argi],"--maxallele") == 0) {
      if(argi+1 == argc || !parse_entire_size(argv[argi+1], &max_allele_len) ||
         max_allele_len == 0)  {
        cmd_print_usage("--maxallele <col> requires an +ve integer argument");
      }
      argi++;
    }
    else if(strcmp(argv[argi],"--maxflank") == 0) {
      if(argi+1 == argc || !parse_entire_size(argv[argi+1], &max_flank_len) ||
         max_flank_len == 0)  {
        cmd_print_usage("--maxflank <col> requires an +ve integer argument");
      }
      argi++;
    }
    else {
      cmd_print_usage("Unknown arg: %s", argv[argi]);
    }
  }

  if(argi + 2 != argc) cmd_print_usage("<out.bubbles.gz> <in.ctx> required");

  char *input_ctx_path = argv[argi];
  char *out_path = argv[argi+1];

  // Open Graph file
  GraphFileReader gfile = INIT_GRAPH_READER;
  int ret = graph_file_open(&gfile, input_ctx_path, false);

  if(ret == 0)
    cmd_print_usage("Cannot read input graph file: %s", input_ctx_path);
  else if(ret < 0)
    cmd_print_usage("Input graph file isn't valid: %s", input_ctx_path);

  // Check reference colours
  for(i = 0; i < num_ref; i++) {
    if(ref_cols[i] >= gfile.hdr.num_of_cols) {
      cmd_print_usage("--ref <col> is greater than max colour [%zu > %u]",
                      ref_cols[i], gfile.hdr.num_of_cols-1);
    }
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

  // Check for compatibility between graph files and path files
  graphs_paths_compatible(&gfile, 1, pfiles, num_pfiles);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t tmp_path_mem, path_mem, thread_mem;
  char path_mem_str[100], thread_mem_str[100];

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +
  // visitedfw/rv(2bits/thread)

  bits_per_kmer = sizeof(Edges)*8 + sizeof(PathIndex)*8 +
                  gfile.hdr.num_of_cols + 2*num_of_threads;

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        gfile.num_of_kmers, false, &graph_mem);

  // Thread memory
  thread_mem = roundup_bits2bytes(kmers_in_hash) * 2;
  bytes_to_str(thread_mem * num_of_threads, 1, thread_mem_str);
  status("[memory] (of which threads: %zu x %zu = %s)\n",
          num_of_threads, thread_mem, thread_mem_str);

  // Path Memory
  tmp_path_mem = path_files_tmp_mem_required(pfiles, num_pfiles);
  path_mem = path_max_mem + tmp_path_mem;

  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  size_t total_mem = graph_mem + thread_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  // Check output file writeable
  if(!futil_is_file_writable(out_path))
    cmd_print_usage("Cannot write output file: %s", out_path);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, gfile.hdr.num_of_cols, 1, kmers_in_hash);

  // Edges merged into one colour
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(uint8_t));

  // In colour
  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);
  db_graph.node_in_cols = calloc2(bytes_per_col*gfile.hdr.num_of_cols, sizeof(uint8_t));

  // Paths
  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(PathIndex));
  memset(db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(PathIndex));

  path_store_alloc(&db_graph.pdata, path_max_mem, tmp_path_mem, path_max_usedcols);

  //
  // Open output file
  //
  gzFile gzout = gzopen(out_path, "w");
  if(gzout == NULL)
    die("Cannot open paths bubble caller output file: %s", out_path);

  //
  // Set up temporary files
  //
  StrBuf *tmppath = strbuf_new();
  char **tmp_paths = malloc2(num_of_threads * sizeof(char*));

  int r = rand() & ((1<<20)-1);

  for(i = 0; i < num_of_threads; i++)
  {
    strbuf_set(tmppath, out_path);
    strbuf_sprintf(tmppath, ".%i.%zu", r, i);
    tmp_paths[i] = strbuf_dup(tmppath);
    if(!futil_is_file_writable(tmp_paths[i])) {
      while(i > 0) unlink(tmp_paths[--i]);
      die("Cannot write temporary file: %s", tmp_paths[i]);
    }
  }
  strbuf_free(tmppath);

  // Load graph
  LoadingStats stats;
  loading_stats_init(&stats);

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  graph_load(&gfile, gprefs, &stats);
  hash_table_print_stats(&db_graph.ht);

  // Load path files
  paths_format_merge(pfiles, num_pfiles, false, &db_graph);

  // Now call variants
  bubble_caller_print_header(&db_graph, gzout, out_path, args);
  invoke_bubble_caller(&db_graph, gzout, num_of_threads, tmp_paths,
                       max_allele_len, max_flank_len, ref_cols, num_ref);

  status("  saved to: %s\n", out_path);
  gzclose(gzout);

  // Clear up threads
  for(i = 0; i < num_of_threads; i++) {
    unlink(tmp_paths[i]);
    free(tmp_paths[i]);
  }
  free(tmp_paths);

  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free(db_graph.kmer_paths);

  path_store_dealloc(&db_graph.pdata);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&gfile);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  return EXIT_SUCCESS;
}
