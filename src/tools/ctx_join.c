
#include "global.h"

#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_format.h"

static const char usage[] =
"usage: ctx_join [options] <mem> <out.ctx> <in1.ctx> [in2.ctx ...]\n"
"  Merge cortex binaries.\n"
"\n"
" Options:\n"
"   --merge    Merge corresponding colours from each binary\n"
"   --flatten  Dump into a single colour binary\n"
"\n"
" Files can be specified with specific colours: samples.ctx:2,3\n";

int main(int argc, char* argv[])
{
  if(argc < 4) print_usage(usage, NULL);

  size_t mem_to_use;
  char *out_ctx_path;
  boolean merge = false, flatten = false;

  int i, argstart;

  for(argstart = 1; argstart < argc; argstart++) {
    if(strcasecmp(argv[argstart],"--merge") == 0) {
      if(merge) warn("merge specified twice");
      merge = true;
    }
    else if(strcasecmp(argv[argstart],"--flatten") == 0) {
      if(flatten) warn("flatten specified twice");
      flatten = true;
    }
    else if(argv[argstart][0] == '-') {
      print_usage(usage, "Unknown argument '%s'", argv[argstart]);
    }
    else break;
  }

  if(argstart + 3 > argc) print_usage(usage, NULL);

  if(!mem_to_integer(argv[argstart], &mem_to_use))
    print_usage(usage, "Invalid <mem> arg (try 1GB or 2M): %s", argv[argstart]);
  argstart++;

  out_ctx_path = argv[argstart++];

  // argstart .. argend-1 are binaries to load
  int num_binaries = argc - argstart;

  // Check all binaries are valid binaries with matching kmer size
  boolean is_binary = false;
  uint32_t kmer_size, kmer_size2, num_of_cols, max_cols, sum_cols;
  uint64_t num_kmers, max_kmers;
  char *ctx_path;

  for(i = 0; i < num_binaries; i++)
  {
    ctx_path = argv[argstart+i];

    if(!binary_probe(ctx_path, &is_binary, &kmer_size2, &num_of_cols, &num_kmers))
      print_usage(usage, "Cannot read input binary file: %s", ctx_path);
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", ctx_path);

    if(i == 0)
    {
      kmer_size = kmer_size2;
      max_kmers = num_kmers;
      max_cols = num_of_cols;
      sum_cols = num_of_cols;
    }
    else
    {
      if(kmer_size != kmer_size2)
        print_usage(usage, "Kmer sizes don't match [%u vs %u]", kmer_size, kmer_size2);

      max_kmers = MAX2(num_kmers, max_kmers);
      max_cols = MAX2(num_of_cols, max_cols);
      sum_cols += num_of_cols;
    }
  }

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  uint32_t cols_used;

  if(flatten) cols_used = 1;
  else if(merge) cols_used = max_cols;
  else cols_used = sum_cols;

  // Pick hash table size
  size_t mem_per_kmer, kmers_in_hash, hash_mem, graph_mem;

  mem_per_kmer = sizeof(BinaryKmer) + (sizeof(Covg) + sizeof(Edges)) * cols_used;
  hash_mem = hash_table_mem2(mem_to_use / mem_per_kmer, &kmers_in_hash);
  graph_mem = kmers_in_hash * mem_per_kmer;

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, cols_used, kmers_in_hash);
  db_graph.col_edges = calloc(db_graph.ht.capacity * cols_used, sizeof(Edges));
  db_graph.col_covgs = calloc(db_graph.ht.capacity * cols_used, sizeof(Covg));

  // Print mem usage
  char graph_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  message("[memory]  graph: %s\n", graph_mem_str);
  hash_table_print_stats(&db_graph.ht);

  message("Using kmer size %u with %u colours\n", kmer_size, cols_used);

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = flatten,
                           .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_colour = -1,
                           .empty_colours = num_binaries == 1,
                           .load_as_union = false,
                           .update_ginfo = false,
                           .db_graph = &db_graph};

  for(i = 0; i < num_binaries; i++)
  {
    if(merge || flatten) prefs.into_colour = 0;
    binary_load(argv[argstart+i], &db_graph, &prefs, stats);
  }

  uint64_t nodes_dumped = binary_dump_graph(out_ctx_path, &db_graph,
                                            CURR_CTX_VERSION, NULL, 0, cols_used);

  hash_table_print_stats(&db_graph.ht);
  printf("Dumped %zu kmers in %u colours\n", (size_t)nodes_dumped, cols_used);
  printf("Done.\n");

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);
}
