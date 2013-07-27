#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_format.h"

static const char usage[] =
"usage: "CMD" intersect [-m <mem>] <graph.ctx> <in.ctx> <out.ctx>\n"
"  Dumps nodes from <in.ctx> that are in <graph.ctx> to <out.ctx>\n";

static void graphs_intersect(const char *graph_ctx_path,
                             const char *in_ctx_path,
                             const char *out_ctx_path,
                             dBGraph *db_graph)
{
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                           .load_seq = false,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_colour = -1,
                           .empty_colours = true,
                           .update_ginfo = false,
                           .db_graph = db_graph};

  binary_load(graph_ctx_path, db_graph, &prefs, stats, NULL);

  // Dump nodes that were flagged
  size_t nodes_dumped = db_graph_filter_file(db_graph, in_ctx_path, out_ctx_path);

  message("Dumped %zu kmers\n", nodes_dumped);
  message("Done.\n");

  seq_loading_stats_free(stats);
}

int ctx_intersect(CmdArgs *args)
{
  cmd_accept_options(args, "m");
  cmd_require_options(args, "m");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 3) print_usage(usage, NULL);

  uint64_t mem_to_use = args->mem_to_use;

  char *graph_ctx_path, *in_ctx_path, *out_ctx_path;

  graph_ctx_path = argv[2];
  in_ctx_path = argv[3];
  out_ctx_path = argv[4];

  // Check graph_ctx_path and in_ctx_path are valid binaries with matching
  // kmer size

  // Probe binary to get kmer_size
  boolean is_binary = false;
  uint32_t kmer_size, kmer_size2, num_of_cols, max_col;
  uint64_t num_kmers, num_kmers2;

  if(!binary_probe(graph_ctx_path, &is_binary, &kmer_size, &num_of_cols, &max_col, &num_kmers))
    print_usage(usage, "Cannot read input binary file: %s", graph_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", graph_ctx_path);

  if(!binary_probe(in_ctx_path, &is_binary, &kmer_size2, &num_of_cols, &max_col, &num_kmers2))
    print_usage(usage, "Cannot read input binary file: %s", in_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", in_ctx_path);

  if(kmer_size != kmer_size2)
    print_usage(usage, "Kmer sizes don't match [%u vs %u]", kmer_size, kmer_size2);

  size_t num_of_hash_kmers;
  size_t req_num_kmers = num_kmers*(1.0/IDEAL_OCCUPANCY);
  size_t hash_mem = hash_table_mem(req_num_kmers, &num_of_hash_kmers);

  size_t total_mem = hash_mem + // hash table
                     num_of_hash_kmers * sizeof(Edges); // edges

  if(total_mem >= mem_to_use)
    die("Not enough memory for the graph");

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  // DEV: Print memory stats

  message("Using kmer size: %u\n", kmer_size);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, num_of_hash_kmers);
  db_graph.edges = calloc(db_graph.ht.capacity, sizeof(Edges));

  graphs_intersect(graph_ctx_path, in_ctx_path, out_ctx_path, &db_graph);

  free(db_graph.edges);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
