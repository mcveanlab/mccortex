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

  // Decide on memory
  size_t hash_kmers, req_num_kmers = gheader.num_of_kmers / IDEAL_OCCUPANCY;
  size_t hash_mem = hash_table_mem(req_num_kmers, &hash_kmers);

  size_t graph_mem = hash_mem +
                     hash_kmers * sizeof(Edges) + // edges
                     hash_kmers * sizeof(uint64_t) + // kmer_paths
                     round_bits_to_bytes(hash_kmers) * gheader.num_of_cols + // in col
                     round_bits_to_bytes(hash_kmers) + // used in contig
                     round_bits_to_bytes(hash_kmers) * 2; // visited fw/rv

  size_t path_mem = input_paths_file == NULL ? 0 : pheader.num_path_bytes;

  // memory to strings
  char graph_mem_str[100], path_mem_str[100], total_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(path_mem, 1, path_mem_str);
  bytes_to_str(graph_mem+path_mem, 1, total_mem_str);

  if(graph_mem + path_mem > args->mem_to_use)
    print_usage(usage, "Not enough memory; requires %s", total_mem_str);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gheader.kmer_size, gheader.num_of_cols, hash_kmers);

  message("[memory]  graph: %s;  paths: %s\n", graph_mem_str, path_mem_str);

  // Edges
  db_graph.edges = calloc2(hash_kmers, sizeof(uint8_t));

  // In colour - used is traversal
  size_t words64_per_col = round_bits_to_words64(hash_kmers);
  db_graph.node_in_cols = calloc2(words64_per_col*gheader.num_of_cols, sizeof(uint64_t));

  // Used in contig
  uint64_t *used_in_contig = calloc2(words64_per_col, sizeof(uint64_t));

  // Paths
  db_graph.kmer_paths = malloc2(hash_kmers * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, hash_kmers * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  path_store_init(&db_graph.pdata, path_store, path_mem, gheader.num_of_cols);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = false,
                           .boolean_covgs = false,
                           .load_seq = false,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_graph = false,
                           .empty_colours = false,
                           .db_graph = &db_graph};

  graph_load(input_ctx_path, &prefs, stats, NULL);
  hash_table_print_stats(&db_graph.ht);

  if(input_paths_file != NULL) {
    paths_format_read(input_paths_file, &pheader, &db_graph,
                      &db_graph.pdata, false);
  }


  // DEV: dump contigs
  // Use kmers unique to the colour with low covg

  free(used_in_contig);
  free(db_graph.edges);
  free(db_graph.node_in_cols);
  free((void *)db_graph.kmer_paths);
  free(path_store);

  graph_header_dealloc(&gheader);
  paths_header_dealloc(&pheader);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  message("  Contigs written to: %s\n", output_fa_path);
  message("Done.\n");

  return EXIT_SUCCESS;
}
