#include "global.h"
#include <time.h>
#include <pthread.h>

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "add_read_paths.h"
#include "binary_format.h"
#include "path_format.h"
#include "graph_walker.h"
#include "bubble_caller.h"

static const char usage[] =
"usage: "CMD" call [-m <mem>|-t <threads>|-p <paths.ctp>] <in.ctx> <out.bubbles.gz>\n"
"  Find bubbles (potential variants) in graph file in.ctx.\n";
// "usage: "CMD" call [options] <out.bubbles.gz> <in.ctx> [in2.ctx ...]\n"
// "  Options:\n"
// "    -m <mem>\n"
// "    -h <kmers>\n";

int ctx_call(CmdArgs *args)
{
  cmd_accept_options(args, "tmp");
  // cmd_require_options(args, "tm", usage);

  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 2) print_usage(usage, NULL);

  uint32_t num_of_threads = args->num_threads;

  char *input_ctx_path = argv[0];
  char *out_path = argv[1];

  // Probe binary to get kmer_size
  boolean is_binary = false;
  GraphFileHeader gheader = {.capacity = 0};

  if(!graph_file_probe(input_ctx_path, &is_binary, &gheader))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  // probe paths file
  const char *input_paths_file = args->ctp_set ? args->ctp_path : NULL;

  boolean valid_paths_file = false;
  PathFileHeader pheader = {.capacity = 0};

  if(input_paths_file != NULL)
  {
    if(!paths_file_probe(input_paths_file, &valid_paths_file, &pheader))
      print_usage(usage, "Cannot read .ctp file: %s", input_paths_file);
    if(!valid_paths_file)
      die("Invalid .ctp file: %s", input_paths_file);
    if(pheader.num_of_cols != gheader.num_of_cols)
      die("Number of colours in .ctp does not match .ctx");
    if(pheader.kmer_size != gheader.kmer_size)
      die("Kmer size in .ctp does not match .ctx");
  }

  // Decide on memory
  size_t kmers_in_hash, ideal_capacity = gheader.num_of_kmers/IDEAL_OCCUPANCY;
  size_t req_num_kmers = args->num_kmers_set ? args->num_kmers : ideal_capacity;
  size_t hash_mem = hash_table_mem(req_num_kmers, &kmers_in_hash);

  size_t graph_mem = hash_mem +
                     kmers_in_hash * sizeof(Edges) + // edges
                     kmers_in_hash * sizeof(uint64_t) + // kmer_paths
                     round_bits_to_bytes(kmers_in_hash) * gheader.num_of_cols + // in col
                     round_bits_to_bytes(kmers_in_hash) * 2; // visited fw/rv

  size_t thread_mem = round_bits_to_bytes(kmers_in_hash) * 2 * num_of_threads;

  // size_t path_mem = mem_to_use - graph_mem - thread_mem;
  size_t path_mem = input_paths_file == NULL ? 0 : pheader.num_path_bytes;

  char graph_mem_str[100], per_thread_mem_str[100], path_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(thread_mem / num_of_threads, 1, per_thread_mem_str);
  bytes_to_str(path_mem, 1, path_mem_str);

  message("[memory]  graph: %s;  threads: %i x %s;  paths: %s\n",
          graph_mem_str, num_of_threads, per_thread_mem_str, path_mem_str);

  if(args->mem_to_use_set && graph_mem+thread_mem > args->mem_to_use)
    die("Not enough memory (increase -m <mem>)");

  if(!test_file_writable(out_path))
    print_usage(usage, "Cannot write output file: %s", out_path);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gheader.kmer_size, gheader.num_of_cols, kmers_in_hash);

  // Edges
  db_graph.edges = calloc2(kmers_in_hash, sizeof(uint8_t));

  // In colour
  size_t words64_per_col = round_bits_to_words64(kmers_in_hash);
  db_graph.node_in_cols = calloc2(words64_per_col*gheader.num_of_cols, sizeof(uint64_t));

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  binary_paths_init(&db_graph.pdata, path_store, path_mem, gheader.num_of_cols);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = false,
                           .boolean_covgs = false,
                           .load_binaries = true,
                           .must_exist_in_graph = false,
                           .empty_colours = true,
                           .db_graph = &db_graph};

  binary_load(input_ctx_path, &prefs, stats, NULL);
  hash_table_print_stats(&db_graph.ht);

  // Load path file
  if(input_paths_file != NULL) {
    paths_format_read(input_paths_file, &pheader, &db_graph,
                      &db_graph.pdata, false);
  }

  /* initialize random seed: */
  srand(time(NULL));

  //
  // Set up temporary files
  //
  StrBuf *tmppath = strbuf_new();
  char **tmp_paths = malloc2(num_of_threads * sizeof(char*));

  int r = rand() & ((1<<20)-1);
  size_t i;

  for(i = 0; i < num_of_threads; i++)
  {
    strbuf_set(tmppath, out_path);
    strbuf_sprintf(tmppath, ".%i.%zu", r, i);
    tmp_paths[i] = strbuf_dup(tmppath);
    if(!test_file_writable(tmp_paths[i])) {
      while(i > 0) unlink(tmp_paths[--i]);
      die("Cannot write temporary file: %s", tmp_paths[i]);
    }
  }

  #ifdef DEBUG
    // db_graph_dump_paths_by_kmer(&db_graph);
  #endif

  // Now call variants
  invoke_bubble_caller(&db_graph, out_path, num_of_threads, tmp_paths);

  // Clear up threads
  for(i = 0; i < num_of_threads; i++) {
    unlink(tmp_paths[i]);
    free(tmp_paths[i]);
  }
  free(tmp_paths);

  free(db_graph.edges);
  free(db_graph.node_in_cols);
  free((void *)db_graph.kmer_paths);
  free(path_store);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  message("Done.\n");
  pthread_exit(NULL);
  return EXIT_SUCCESS;
}
