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
  size_t num_binaries = argc - argi;
  char **binary_paths = argv + argi;

  // Check all binaries are valid binaries with matching kmer size
  boolean is_binary = false;
  size_t i, col, kmer_size = 0, ctx_max_kmers = 0;
  uint32_t ctx_max_cols[num_binaries], ctx_num_cols[num_binaries];
  GraphFileHeader gheader = {.capacity = 0};
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

    ctx_num_cols[i] = gheader.num_of_cols;
    ctx_max_cols[i] = gheader.max_col;
    ctx_max_kmers = MAX2(gheader.num_of_kmers, ctx_max_kmers);

    if(i == 0)
      kmer_size = gheader.kmer_size;
    else if(kmer_size != gheader.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%zu vs %u]",
                  kmer_size, gheader.kmer_size);
    }
  }

  // Probe intersection files
  uint64_t min_intersect_num_kmers = 0;

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

  if(num_intersect > 0)
    ctx_max_kmers = min_intersect_num_kmers;

  //
  // Decide on memory
  //
  size_t extra_bits_per_kmer, kmers_in_hash;

  extra_bits_per_kmer = (sizeof(Covg) + sizeof(Edges)) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        ctx_max_kmers, true);

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity, sizeof(Covg));

  // Load intersection binaries
  char *intsct_gname_ptr = NULL;
  StrBuf intersect_gname;
  strbuf_alloc(&intersect_gname, 1024);

  if(num_intersect > 0)
  {
    intsct_gname_ptr = intersect_gname.buff;
    SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                             .boolean_covgs = true, // covg++
                             .load_seq = true,
                             .quality_cutoff = 0, .ascii_fq_offset = 0,
                             .homopolymer_cutoff = 0,
                             .remove_dups_se = false, .remove_dups_pe = false,
                             .load_binaries = true,
                             .must_exist_in_graph = false,
                             .empty_colours = false,
                             .db_graph = &db_graph};

    for(i = 0; i < num_intersect; i++)
    {
      prefs.must_exist_in_graph = (i > 0);

      graph_load(intersect_paths[i], &prefs, NULL, &gheader);

      // Update intersect header
      for(col = 0; col < gheader.num_of_cols; col++)
        graph_info_make_intersect(&gheader.ginfo[col], &intersect_gname);
    }

    if(num_intersect > 1)
    {
      // Remove nodes where covg != num_intersect
      HASH_TRAVERSE(&db_graph.ht, remove_non_intersect_nodes, db_graph.col_covgs,
                    num_intersect, &db_graph.ht);
    }

    // Zero covgs and edges
    memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
    memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));

    status("Loaded intersection set\n");
  }

  graph_files_merge(out_ctx_path, binary_paths, num_binaries,
                    ctx_num_cols, ctx_max_cols,
                    merge, flatten, intsct_gname_ptr, &db_graph);

  graph_header_dealloc(&gheader);
  strbuf_dealloc(&intersect_gname);

  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
