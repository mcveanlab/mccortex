#include "global.h"
#include <time.h>

#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_paths.h"
#include "binary_format.h"
#include "path_format.h"

static const char usage[] =
"usage: ctx_contigs <mem> <in.ctx> <colour> <out.fa>\n"
"  Pull out contigs for a given colour\n";

int main(int argc, char* argv[])
{
  if(argc != 5) print_usage(usage, NULL);

  char *input_ctx_path = argv[2], *output_fa_path = argv[4];

  size_t mem_to_use = 0;
  uint32_t colour;

  if(!mem_to_integer(argv[1], &mem_to_use) || mem_to_use == 0)
    print_usage(usage, "Invalid memory argument: %s", argv[1]);

  // Probe ctx
  boolean is_binary = false;
  uint32_t kmer_size, num_of_cols;
  uint64_t num_kmers;

  if(!binary_probe(input_ctx_path, &is_binary, &kmer_size, &num_of_cols, &num_kmers))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  if(!parse_entire_uint(argv[3], &colour) || colour == 0)
    print_usage(usage, "Invalid colour: %s", argv[3]);

  if(!test_file_writable(output_fa_path))
    print_usage(usage, "Cannot write to output file: %s", output_fa_path);

  // DEV: look for ctp file
  size_t path_mem = 0;

  // Decide on memory
  size_t hash_kmers, req_num_kmers = num_kmers*(1.0/IDEAL_OCCUPANCY);
  size_t hash_mem = hash_table_mem(req_num_kmers, &hash_kmers);

  size_t graph_mem = hash_mem +
                     hash_kmers * sizeof(Edges) + // edges
                     hash_kmers * sizeof(uint64_t) * 2 + // kmer_paths
                     round_bits_to_bytes(hash_kmers) + // in col
                     round_bits_to_bytes(hash_kmers) * 2; // visited fw/rv

  // memory to strings
  char graph_mem_str[100], path_mem_str[100], total_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(path_mem, 1, path_mem_str);
  bytes_to_str(graph_mem+path_mem, 1, total_mem_str);

  if(graph_mem > mem_to_use)
    print_usage(usage, "Not enough memory; requires %s", total_mem_str);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, num_of_cols, hash_kmers);

  message("[memory]  graph: %s;  paths: %s\n", graph_mem_str, path_mem_str);

  // Edges
  db_graph.edges = calloc(hash_kmers, sizeof(uint8_t));
  if(db_graph.edges == NULL) die("Out of memory");

  // In colour - used is traversal
  size_t words64_per_col = round_bits_to_words64(hash_kmers);
  db_graph.node_in_cols = calloc(words64_per_col*num_of_cols, sizeof(uint64_t));
  if(db_graph.node_in_cols == NULL) die("Out of memory");

  // Paths
  db_graph.kmer_paths = malloc(hash_kmers * sizeof(uint64_t) * 2);
  if(db_graph.kmer_paths == NULL) die("Out of memory");
  memset(db_graph.kmer_paths, 0xff, hash_kmers * sizeof(uint64_t) * 2);

  uint8_t *path_store = malloc(path_mem);
  if(path_store == NULL) die("Out of memory");
  binary_paths_init(&db_graph.pdata, path_store, path_mem, num_of_cols);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = false,
                           .load_seq = false,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_colour = -1,
                           .empty_colours = false,
                           .update_ginfo = true,
                           .db_graph = &db_graph};

  binary_load(input_ctx_path, &db_graph, &prefs, stats);
  hash_table_print_stats(&db_graph.ht);


  // DEV: dump contigs

  free(db_graph.edges);
  free(db_graph.node_in_cols);
  free(db_graph.kmer_paths);
  free(path_store);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  message("  Contigs written to: %s\n", output_fa_path);
  message("Done.\n");
}
