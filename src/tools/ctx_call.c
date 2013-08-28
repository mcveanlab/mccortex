#include "global.h"
#include <time.h>
#include <pthread.h>
#include <unistd.h>

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "add_read_paths.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_walker.h"
#include "bubble_caller.h"

// "usage: "CMD" call [options] <out.bubbles.gz> <in.ctx> [in2.ctx ...]\n"
static const char usage[] =
"usage: "CMD" call [-m <mem>|-t <threads>|-p <paths.ctp>] <in.ctx> <out.bubbles.gz>\n"
"  Find bubbles (potential variants) in graph file in.ctx.\n"
"  Options:  -m <mem> | -h <kmers> | -t <threads> | -p <paths.ctp>\n";

int ctx_call(CmdArgs *args)
{
  cmd_accept_options(args, "tmp");

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
  if(args->num_ctp_files > 1)
    print_usage(usage, "Sorry, we only accept one -p <in.ctp> at the moment");

  const char *input_paths_file = args->num_ctp_files ? args->ctp_files[0] : NULL;

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

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, thread_mem, total_mem;
  char path_mem_str[100], thread_mem_str[100], total_mem_str[100];

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +
  // visitedfw/rv(2bits/kmer/thread)

  bits_per_kmer = sizeof(Edges)*8 + sizeof(uint64_t)*8 +
                  gheader.num_of_cols + 2*num_of_threads;

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        gheader.num_of_kmers, true);

  // Thread memory
  thread_mem = round_bits_to_bytes(kmers_in_hash) * 2;
  bytes_to_str(thread_mem * num_of_threads, 1, thread_mem_str);
  status("[memory] (of which threads: %u x %zu = %s)\n",
          num_of_threads, thread_mem, thread_mem_str);

  // Path Memory
  bytes_to_str(pheader.num_path_bytes, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  total_mem = pheader.num_path_bytes + (kmers_in_hash*bits_per_kmer)/8;
  bytes_to_str(total_mem, 1, total_mem_str);

  if(total_mem > args->mem_to_use)
    die("Requires at least %s memory", total_mem_str);

  // Check output file writeable
  if(!test_file_writable(out_path))
    print_usage(usage, "Cannot write output file: %s", out_path);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gheader.kmer_size, gheader.num_of_cols, kmers_in_hash);

  if(kmers_in_hash != db_graph.ht.capacity) die("Mismatch");

  // Edges
  db_graph.edges = calloc2(kmers_in_hash, sizeof(uint8_t));

  // In colour
  size_t words64_per_col = round_bits_to_words64(kmers_in_hash);
  db_graph.node_in_cols = calloc2(words64_per_col*gheader.num_of_cols, sizeof(uint64_t));

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  uint8_t *path_store = malloc2(pheader.num_path_bytes);
  path_store_init(&db_graph.pdata, path_store, pheader.num_path_bytes, gheader.num_of_cols);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = false,
                           .boolean_covgs = false,
                           .load_binaries = true,
                           .must_exist_in_graph = false,
                           .empty_colours = true,
                           .db_graph = &db_graph};

  graph_load(input_ctx_path, &prefs, stats, NULL);
  hash_table_print_stats(&db_graph.ht);

  // Load path file
  if(input_paths_file != NULL) {
    paths_format_read(input_paths_file, &pheader, &db_graph,
                      &db_graph.pdata, false);
  }

  /* initialize random seed: */
  srand(time(NULL) + getpid());

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
  strbuf_free(tmppath);

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

  graph_header_dealloc(&gheader);
  paths_header_dealloc(&pheader);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
