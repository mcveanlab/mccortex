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
#include "graph_file_filter.h"

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
"   -m <mem>      Memory to use\n"
"   -h <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   --usecols <c> How many colours to load at once [default: 1]\n"
"   --merge       Merge corresponding colours from each graph file\n"
"   --flatten     Dump into a single colour graph\n"
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
  size_t num_intersect = 0, use_ncols = 1;

  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++) {
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
        print_usage(usage, "--intersect <A.ctx> requires an argument");
      num_intersect++;
      argi++;
    }
    else if(strcasecmp(argv[argi],"--usecols") == 0)
    {
      unsigned int tmp;
      if(argi+1 >= argc || !parse_entire_uint(argv[argi+1], &tmp))
        print_usage(usage, "--usecols <c> requires a positive integer argument");
      use_ncols = tmp;
      argi++;
    }
    else {
      print_usage(usage, "Unknown argument '%s'", argv[argi]);
    }
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
  size_t num_graphs = argc - argi;
  char **paths = argv + argi;

  status("Probing %zu graph files", num_graphs);

  // Check all binaries are valid binaries with matching kmer size
  size_t i, col, ncols, ctx_max_kmers = 0, max_cols = 0, sum_cols = 0, total_cols;
  size_t min_intersect_num_kmers = 0;
  GraphFileReader files[num_graphs], intersect_files[num_intersect];

  for(i = 0; i < num_graphs; i++)
  {
    files[i] = INIT_GRAPH_READER;
    int ret = graph_file_open(&files[i], paths[i], false);

    if(ret == 0)
      print_usage(usage, "Cannot read input binary file: %s", paths[i]);
    else if(ret < 0)
      print_usage(usage, "Input binary file isn't valid: %s", paths[i]);

    if(i > 0 && files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, files[i].hdr.kmer_size);
    }

    if(flatten) {
      files[i].flatten = true;
      files[i].intocol = 0;
    }

    ncols = files[i].intocol + graph_file_outncols(&files[i]);
    max_cols = MAX2(max_cols, ncols);
    sum_cols += ncols;
    ctx_max_kmers = MAX2(ctx_max_kmers, files[i].hdr.num_of_kmers);
  }

  if(flatten) total_cols = 1;
  else if(merge) total_cols = max_cols;
  else {
    total_cols = 0;
    for(i = 0; i < num_graphs; i++) {
      size_t offset = total_cols;
      total_cols += files[i].intocol + graph_file_outncols(&files[i]);
      files[i].intocol += offset;
    }
  }

  // Probe intersection files
  for(i = 0; i < num_intersect; i++)
  {
    intersect_files[i] = INIT_GRAPH_READER;
    int ret = graph_file_open(&intersect_files[i], intersect_paths[i], false);

    if(ret == 0)
      print_usage(usage, "Cannot read input binary file: %s", intersect_paths[i]);
    else if(ret < 0)
      print_usage(usage, "Input binary file isn't valid: %s", intersect_paths[i]);

    if(files[0].hdr.kmer_size != intersect_files[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, intersect_files[i].hdr.kmer_size);
    }

    intersect_files[i].flatten = true;

    if(i == 0) min_intersect_num_kmers = intersect_files[i].hdr.num_of_kmers;
    else if(intersect_files[i].hdr.num_of_kmers < min_intersect_num_kmers)
    {
      // Put smallest intersection binary first
      GraphFileReader tmp;
      SWAP(intersect_files[i], intersect_files[0], tmp);
      min_intersect_num_kmers = intersect_files[i].hdr.num_of_kmers;
    }
  }

  if(num_intersect > 0)
    ctx_max_kmers = min_intersect_num_kmers;

  if(use_ncols > 1 && flatten) {
    warn("I only need one colour for '--flatten' ('--usecols %zu' ignored)", use_ncols);
    use_ncols = 1;
  }
  else if(use_ncols > total_cols) {
    warn("I only need %zu colour%s ('--usecols %zu' ignored)",
         total_cols, (total_cols != 1 ? "s" : ""), use_ncols);
    use_ncols = total_cols;
  }

  status("Output %zu cols; from %zu files; intersecting %zu graphs; using %zu cols in memory",
         total_cols, num_graphs, num_intersect, use_ncols);

  if(num_graphs == 1 && num_intersect == 0)
  {
    // Loading only one binary with no intersection filter
    // don't need to store a graph in memory
    graph_stream_filter2(out_ctx_path, &files[0], NULL, NULL);
    graph_file_dealloc(&files[0]);
    return EXIT_SUCCESS;
  }

  //
  // Decide on memory
  //
  size_t extra_bits_per_kmer, kmers_in_hash;

  extra_bits_per_kmer = (sizeof(Covg) + sizeof(Edges)) * 8 * use_ncols;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        ctx_max_kmers, true);

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, use_ncols, use_ncols, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity*use_ncols, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity*use_ncols, sizeof(Covg));

  // Load intersection binaries
  char *intsct_gname_ptr = NULL;
  StrBuf intersect_gname;
  strbuf_alloc(&intersect_gname, 1024);

  if(num_intersect > 0)
  {
    intsct_gname_ptr = intersect_gname.buff;
    SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                             // binaries
                             .boolean_covgs = true, // covg++ only
                             .must_exist_in_graph = false,
                             .empty_colours = false,
                             // sequence
                             .quality_cutoff = 0, .ascii_fq_offset = 0,
                             .homopolymer_cutoff = 0,
                             .remove_dups_se = false, .remove_dups_pe = false};

    Edges *tmpcoledges = db_graph.col_edges;
    db_graph.col_edges = NULL;

    for(i = 0; i < num_intersect; i++)
    {
      prefs.must_exist_in_graph = (i > 0);

      graph_load2(&intersect_files[i], &prefs, NULL);

      // Update intersect header
      ncols = graph_file_outncols(&intersect_files[i]);
      for(col = 0; col < ncols; col++)
        graph_info_make_intersect(&intersect_files[i].hdr.ginfo[col], &intersect_gname);
    }

    if(num_intersect > 1)
    {
      // Remove nodes where covg != num_intersect
      HASH_TRAVERSE(&db_graph.ht, remove_non_intersect_nodes, db_graph.col_covgs,
                    num_intersect, &db_graph.ht);
    }

    // Zero covgs
    memset(db_graph.col_covgs, 0, db_graph.ht.capacity * use_ncols * sizeof(Covg));

    // Swap edges with NULL instead of having to wipe
    db_graph.col_edges = tmpcoledges;

    for(i = 0; i < num_intersect; i++) graph_file_dealloc(&intersect_files[i]);

    // Reset graph info
    for(i = 0; i < db_graph.num_of_cols; i++)
      graph_info_init(&db_graph.ginfo[i]);

    status("Loaded intersection set\n");
  }

  graph_files_merge2(out_ctx_path, files, num_graphs,
                     intsct_gname_ptr, &db_graph);

  for(i = 0; i < num_graphs; i++) graph_file_dealloc(&files[i]);

  strbuf_dealloc(&intersect_gname);

  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
