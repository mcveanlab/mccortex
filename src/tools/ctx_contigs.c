#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "path_store.h"
#include "graph_format.h"
#include "path_format.h"

static const char usage[] =
"usage: "CMD" contigs [-m <mem>|-h <kmers>|-p <paths>] <in.ctx> <colour> <out.fa>\n"
"  Pull out contigs for a given colour\n";

int ctx_contigs(CmdArgs *args)
{
  cmd_accept_options(args, "mhp");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  char *input_ctx_path = argv[2], *output_fa_path = argv[4];
  uint32_t colour;

  // Probe ctx
  boolean is_binary = false;
  GraphFileHeader gheader = {.capacity = 0};

  if(!graph_file_probe(input_ctx_path, &is_binary, &gheader))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  if(!parse_entire_uint(argv[3], &colour) || colour == 0)
    print_usage(usage, "Invalid colour: %s", argv[3]);

  if(!test_file_writable(output_fa_path))
    print_usage(usage, "Cannot write to output file: %s", output_fa_path);

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
  size_t bits_per_kmer, kmers_in_hash, total_mem;
  char path_mem_str[100], total_mem_str[100];

  // edges, kmer_paths, in_colour, visited(fw/rv) [2bits], used in contig [1bit]
  bits_per_kmer = sizeof(Edges)*8 + sizeof(uint64_t)*8 + gheader.num_of_cols + 3;

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        gheader.num_of_kmers, true);

  // Path Memory
  bytes_to_str(pheader.num_path_bytes, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  total_mem = pheader.num_path_bytes + (kmers_in_hash*bits_per_kmer)/8;
  bytes_to_str(total_mem, 1, total_mem_str);

  if(total_mem > args->mem_to_use)
    die("Requires at least %s memory", total_mem_str);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gheader.kmer_size, gheader.num_of_cols, 1, kmers_in_hash);

  // Edges
  db_graph.col_edges = calloc2(kmers_in_hash, sizeof(uint8_t));

  // In colour - used is traversal
  size_t words64_per_col = round_bits_to_words64(kmers_in_hash);
  db_graph.node_in_cols = calloc2(words64_per_col*gheader.num_of_cols, sizeof(uint64_t));

  // Used in contig
  uint64_t *used_in_contig = calloc2(words64_per_col, sizeof(uint64_t));

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  uint8_t *path_store = malloc2(pheader.num_path_bytes);
  path_store_init(&db_graph.pdata, path_store,
                  pheader.num_path_bytes, gheader.num_of_cols);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .db_graph = &db_graph,
                           .merge_colours = false,
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = false};

  graph_load(input_ctx_path, &prefs, stats, NULL);
  hash_table_print_stats(&db_graph.ht);

  if(input_paths_file != NULL) {
    paths_format_read(input_paths_file, &pheader, &db_graph,
                      &db_graph.pdata, false);
  }


  // DEV: dump contigs
  // Use kmers unique to the colour with low covg

  free(used_in_contig);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free((void *)db_graph.kmer_paths);
  free(path_store);

  graph_header_dealloc(&gheader);
  paths_header_dealloc(&pheader);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  status("Contigs written to: %s\n", output_fa_path);

  return EXIT_SUCCESS;
}
