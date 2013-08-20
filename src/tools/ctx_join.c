#include "global.h"

#include "string_buffer.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "graph_info.h"
#include "db_node.h"
#include "graph_format.h"

// Given (A,B,C) are ctx binaries, A:1 means colour 1 in A,
// {A:1,B:0} is loading A:1 and B:0 into a single colour
//
// Default behaviour is to load colours consecutively
//   output: {A:0},{A:1},{B:0},{B:1},{B:2},{C:0},{C:1}
//
// --flatten
//   All colours into one
//   output: {A:0,A:1,B:0,B:1,B:2,C:0,C:1}
//
// --merge
//   Join input colours
//   output: {A:0,B:0,C:0},{A:1,B:1,C:1},{B:2}

static const char usage[] =
"usage: "CMD" join [options] <out.ctx> [offset:]in1.ctx[:1,2,4-5] [in2.ctx ...]\n"
"  Merge cortex binaries.  \n"
"\n"
"  Options:\n"
"   -m <mem>     Memory to use\n"
"   -h <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   --merge      Merge corresponding colours from each binary\n"
"   --flatten    Dump into a single colour binary\n"
"\n"
"   --intersect <a.ctx>\n"
"     Only load the kmers that are in graph A.ctx. Can be specified multiple times.\n"
"     <a.ctx> is NOT merged into the output file.\n"
"\n"
"  Files can be specified with specific colours: samples.ctx:2,3\n"
"  Offset specifies where to load the first colour (merge only).\n";

static inline void remove_non_intersect_nodes(hkey_t node, Covg *covgs,
                                              Covg num, HashTable *ht)
{
  if(covgs[node] != num)
    hash_table_delete(ht, node);
}

int ctx_join(CmdArgs *args)
{
  cmd_accept_options(args, "mh");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  char *out_ctx_path;
  boolean merge = false, flatten = false;

  int argi;
  size_t num_intersect = 0;

  for(argi = 0; argi < argc; argi++) {
    if(strcasecmp(argv[argi],"--merge") == 0) {
      if(merge) warn("merge specified twice");
      merge = true;
    }
    else if(strcasecmp(argv[argi],"--flatten") == 0) {
      if(flatten) warn("flatten specified twice");
      flatten = true;
    }
    else if(strcasecmp(argv[argi],"--intersect") == 0)
    {
      if(argi+1 >= argc)
        print_usage(usage, "--intersect <A.ctx> requires and argument");
      num_intersect++;
      argi++;
    }
    else if(argv[argi][0] == '-') {
      print_usage(usage, "Unknown argument '%s'", argv[argi]);
    }
    else break;
  }

  // Store intersection files
  char *intersect_paths[num_intersect+1];
  int argj;
  num_intersect = 0;

  for(argj = 0; argj < argi; argj++) {
    if(strcasecmp(argv[argj],"--intersect") == 0) {
      intersect_paths[num_intersect++] = argv[argj+1];
      argj++;
    }
  }

  if(argc - argi < 2)
    print_usage(usage, "Please specify output and input binaries");

  out_ctx_path = argv[argi++];

  // argi .. argend-1 are binaries to load
  uint32_t num_binaries = argc - argi;
  char **binary_paths = argv + argi;

  // Check all binaries are valid binaries with matching kmer size
  boolean is_binary = false;
  uint32_t i, kmer_size = 0;
  uint32_t ctx_max_cols[num_binaries], ctx_num_cols[num_binaries];
  GraphFileHeader gheader;
  char *path;

  for(i = 0; i < num_binaries; i++)
  {
    // Strip off offset (e.g. 12:in.ctx)
    for(path = binary_paths[i]; *path >= '0' && *path <= '9'; path++) {}

    if(path > binary_paths[i] && *path == ':') path++;
    else path = binary_paths[i];

    if(!graph_file_probe(path, &is_binary, &gheader))
      print_usage(usage, "Cannot read input binary file: %s", binary_paths[i]);
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", binary_paths[i]);

    if(i == 0)
      kmer_size = gheader.kmer_size;
    else if(kmer_size != gheader.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  kmer_size, gheader.kmer_size);
    }
  }

  // Probe intersection files
  uint64_t min_intersect_num_kmers;

  for(i = 0; i < num_intersect; i++)
  {
    if(!graph_file_probe(intersect_paths[i], &is_binary, &gheader))
      print_usage(usage, "Cannot read intersect binary file: %s", intersect_paths[i]);
    else if(!is_binary)
      print_usage(usage, "Intersect binary file isn't valid: %s", intersect_paths[i]);

    if(i == 0) min_intersect_num_kmers = gheader.num_of_kmers;
    else if(gheader.num_of_kmers < min_intersect_num_kmers)
    {
      // Put smallest intersection binary first
      char *tmpstr;
      SWAP(intersect_paths[i], intersect_paths[0], tmpstr);
      min_intersect_num_kmers = gheader.num_of_kmers;
    }
  }

  // Pick hash table size
  size_t kmers_in_hash;
  size_t extra_mem_per_kmer = sizeof(Covg) + sizeof(Edges);
  size_t mem_per_kmer = sizeof(BinaryKmer) + extra_mem_per_kmer;

  if(num_intersect > 0)
  {
    size_t ideal_occupancy = min_intersect_num_kmers / IDEAL_OCCUPANCY;
    if(args->num_kmers_set) {
      char tmp[50];
      if(args->num_kmers < min_intersect_num_kmers) {
        die("Need at least -h %s", bytes_to_str(min_intersect_num_kmers,1,tmp));
      } else if(args->num_kmers < ideal_occupancy) {
        warn("Less than ideal hash occupancy (-h %zu < %zu [70%%])",
             args->num_kmers, ideal_occupancy);
      }
      kmers_in_hash = args->num_kmers;
    }
    else kmers_in_hash = ideal_occupancy;

    size_t hash_mem = hash_table_mem(kmers_in_hash, &kmers_in_hash);

    size_t mem_needed = hash_mem + kmers_in_hash * mem_per_kmer;
    char mem_needed_str[50], mem_set_str[50];
    bytes_to_str(mem_needed, 1, mem_needed_str);
    bytes_to_str(args->mem_to_use, 1, mem_set_str);

    if(args->mem_to_use_set && mem_needed > args->mem_to_use)
      die("Requires %s memory, you set -m %s", mem_needed_str, mem_set_str);

    message("[memory] graph: %s\n", mem_needed_str);
  }
  else
    kmers_in_hash = cmd_get_kmers_in_hash(args, extra_mem_per_kmer);

  message("\n");

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity, sizeof(Covg));

  // Load intersection binaries
  if(num_intersect > 0)
  {
    SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                             .boolean_covgs = true,
                             .load_seq = true,
                             .quality_cutoff = 0, .ascii_fq_offset = 0,
                             .homopolymer_cutoff = 0,
                             .remove_dups_se = false, .remove_dups_pe = false,
                             .load_binaries = true,
                             .must_exist_in_graph = false,
                             .empty_colours = false,
                             .db_graph = &db_graph};

    graph_load(intersect_paths[0], &prefs, NULL, NULL);

    if(num_intersect > 1)
    {
      prefs.must_exist_in_graph = true;

      for(i = 1; i < num_intersect; i++)
        graph_load(intersect_paths[i], &prefs, NULL, NULL);

      // Remove nodes where covg != num_intersect
      HASH_TRAVERSE(&db_graph.ht, remove_non_intersect_nodes, db_graph.col_covgs,
                    num_intersect, &db_graph.ht);
    }

    // Zero covgs and edges
    memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
    memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));

    message("Loaded intersection set\n\n");
  }

  graph_files_merge(out_ctx_path, binary_paths, num_binaries,
                    ctx_num_cols, ctx_max_cols,
                    merge, flatten, num_intersect > 0, &db_graph);

  graph_header_dealloc(&gheader);

  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);
  message("Done.\n");
  return EXIT_SUCCESS;
}
