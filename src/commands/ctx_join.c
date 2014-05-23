#include "global.h"

#include "commands.h"
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
// --overlap
//   All colour 0s go into colour 0, colour 1s go into colour 1 etc.
//   output: {A:0,B:0,C:0},{A:1,B:1,C:1},{B:2}

const char join_usage[] =
"usage: "CMD" join [options] in1.ctx [[offset:]in2.ctx[:1,2,4-5] ...]\n"
"  Merge cortex graphs.\n"
"\n"
"  Options:\n"
"   -o, --out <out.ctp>     Output file [required]\n"
"   -m, --memory <mem>      Memory to use\n"
"   -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   -s, --ncols <c>         How many colours to load at once [default: 1]\n"
"   -v, --overlap           Merge corresponding colours from each graph file\n"
"   -f, --flatten           Dump into a single colour graph\n"
"   -i, --intersect <a.ctx> Only load the kmers that are in graph A.ctx. Can be\n"
"                           specified multiple times. <a.ctx> is NOT merged into\n"
"                           the output file.\n"
"\n"
"  Files can be specified with specific colours: samples.ctx:2,3\n"
"  Offset specifies where to load the first colour: 3:samples.ctx\n";

static inline void remove_non_intersect_nodes(hkey_t node, Covg *covgs,
                                              Covg num, HashTable *ht)
{
  if(covgs[node] != num)
    hash_table_delete(ht, node);
}

int ctx_join(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked we have at least 1 argument

  char *out_ctx_path = args->output_file;
  bool overlap = false, flatten = false;
  char **intersect_paths = ctx_calloc(argc / 2, sizeof(char*));
  GraphFileReader *gfiles = ctx_calloc(argc / 2, sizeof(GraphFileReader));
  GraphFileReader *intersect_gfiles = ctx_calloc(argc / 2, sizeof(GraphFileReader));
  size_t num_intersect = 0;

  int argi;
  size_t use_ncols = args->use_ncols; // may use fewer colours if some not needed

  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++) {
    if(strcmp(argv[argi],"--overlap") == 0 || strcmp(argv[argi],"-v") == 0) {
      if(overlap) warn("overlap specified twice");
      overlap = true;
    }
    else if(strcmp(argv[argi],"--flatten") == 0 || strcmp(argv[argi],"-f") == 0) {
      if(flatten) warn("flatten specified twice");
      flatten = true;
    }
    else if(strcmp(argv[argi],"--intersect") == 0 || strcmp(argv[argi],"-i") == 0)
    {
      if(argi+1 >= argc)
        cmd_print_usage("--intersect <A.ctx> requires an argument");
      intersect_paths[num_intersect++] = argv[argi+1];
      argi++;
    }
    else {
      cmd_print_usage("Unknown argument '%s'", argv[argi]);
    }
  }

  if(argi >= argc)
    cmd_print_usage("Please specify at least one input graph file");

  // argi .. argend-1 are graphs to load
  size_t num_gfiles = (size_t)(argc - argi);
  char **paths = argv + argi;

  status("Probing %zu graph files and %zu intersect files", num_gfiles, num_intersect);

  // Check all binaries are valid binaries with matching kmer size
  size_t i, col, ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;
  size_t max_cols = 0, sum_cols = 0, total_cols;
  size_t min_intersect_num_kmers = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    gfiles[i] = INIT_GRAPH_READER;
    graph_file_open(&gfiles[i], paths[i], true);

    if(gfiles[0].hdr.kmer_size != gfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                  gfiles[0].hdr.kmer_size, gfiles[i].hdr.kmer_size);
    }

    if(flatten) {
      gfiles[i].fltr.flatten = true;
      file_filter_update_intocol(&gfiles[i].fltr, 0);
    }

    ncols = graph_file_usedcols(&gfiles[i]);
    max_cols = MAX2(max_cols, ncols);
    sum_cols += ncols;
    ctx_max_kmers = MAX2(ctx_max_kmers, gfiles[i].num_of_kmers);
    ctx_sum_kmers += gfiles[i].num_of_kmers;
  }

  if(flatten) total_cols = 1;
  else if(overlap) total_cols = max_cols;
  else {
    total_cols = 0;
    for(i = 0; i < num_gfiles; i++) {
      size_t offset = total_cols;
      total_cols += graph_file_usedcols(&gfiles[i]);
      file_filter_update_intocol(&gfiles[i].fltr, gfiles[i].fltr.intocol + offset);
    }
  }

  // Probe intersection graph files
  for(i = 0; i < num_intersect; i++)
  {
    intersect_gfiles[i] = INIT_GRAPH_READER;
    graph_file_open(&intersect_gfiles[i], intersect_paths[i], true);

    if(gfiles[0].hdr.kmer_size != intersect_gfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                  gfiles[0].hdr.kmer_size, intersect_gfiles[i].hdr.kmer_size);
    }

    intersect_gfiles[i].fltr.flatten = true;

    if(i == 0) min_intersect_num_kmers = intersect_gfiles[i].num_of_kmers;
    else if(intersect_gfiles[i].num_of_kmers < min_intersect_num_kmers)
    {
      // Put smallest intersection binary first
      SWAP(intersect_gfiles[i], intersect_gfiles[0]);
      min_intersect_num_kmers = intersect_gfiles[i].num_of_kmers;
    }
  }

  bool take_intersect = (num_intersect > 0);

  if(take_intersect)
    ctx_max_kmers = min_intersect_num_kmers;

  if(use_ncols < total_cols && strcmp(out_ctx_path,"-") == 0) {
    if(args->use_ncols_set)
      die("I need %zu colours if outputting to STDOUT (--ncols)", total_cols);
    else
      use_ncols = total_cols;
  }

  if(use_ncols > 1 && flatten) {
    warn("I only need one colour for '--flatten' ('--ncols %zu' ignored)", use_ncols);
    use_ncols = 1;
  }
  else if(use_ncols > total_cols) {
    warn("I only need %zu colour%s ('--ncols %zu' ignored)",
         total_cols, util_plural_str(total_cols), use_ncols);
    use_ncols = total_cols;
  }

  status("Output %zu cols; from %zu files; intersecting %zu graphs; "
         "using %zu cols in memory",
         total_cols, num_gfiles, num_intersect, use_ncols);

  if(num_gfiles == 1 && num_intersect == 0)
  {
    // Loading only one binary with no intersection filter
    // don't need to store a graph in memory
    graph_stream_filter_mkhdr(out_ctx_path, &gfiles[0], NULL, NULL, NULL);
    graph_file_close(&gfiles[0]);
    return EXIT_SUCCESS;
  }

  //
  // Decide on memory
  //
  size_t extra_bits_per_kmer, kmers_in_hash, graph_mem;

  extra_bits_per_kmer = (sizeof(Covg) + sizeof(Edges)) * 8 * use_ncols;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  // Check out_ctx_path is writable
  if(strcmp(out_ctx_path,"-") != 0 && !futil_is_file_writable(out_ctx_path))
    cmd_print_usage("Cannot write to output: %s", out_ctx_path);

  // Create db_graph
  dBGraph db_graph;
  Edges *intersect_edges = NULL;

  // Edges
  if(take_intersect) {
    db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, 1, 1, kmers_in_hash);
    db_graph.col_edges = ctx_calloc(db_graph.ht.capacity*(use_ncols+1), sizeof(Edges));
  }
  else {
    db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, use_ncols, use_ncols, kmers_in_hash);
    db_graph.col_edges = ctx_calloc(db_graph.ht.capacity*use_ncols, sizeof(Edges));
  }

  db_graph.col_covgs = ctx_calloc(db_graph.ht.capacity*use_ncols, sizeof(Covg));

  // Load intersection binaries
  char *intsct_gname_ptr = NULL;
  StrBuf intersect_gname;
  strbuf_alloc(&intersect_gname, 1024);

  if(take_intersect)
  {
    GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                                .boolean_covgs = true, // covg++ only
                                .must_exist_in_graph = false,
                                .must_exist_in_edges = NULL,
                                .empty_colours = false};

    for(i = 0; i < num_intersect; i++)
    {
      graph_load(&intersect_gfiles[i], gprefs, NULL);

      // Update intersect header
      ncols = graph_file_outncols(&intersect_gfiles[i]);
      for(col = 0; col < ncols; col++) {
        graph_info_make_intersect(&intersect_gfiles[i].hdr.ginfo[col],
                                  &intersect_gname);
      }

      gprefs.must_exist_in_graph = true;
      gprefs.must_exist_in_edges = db_graph.col_edges;
    }

    if(num_intersect > 1)
    {
      // Remove nodes where covg != num_intersect
      HASH_ITERATE_SAFE(&db_graph.ht, remove_non_intersect_nodes,
                        db_graph.col_covgs, (Covg)num_intersect, &db_graph.ht);
    }

    status("Loaded intersection set\n");
    intsct_gname_ptr = intersect_gname.buff;

    for(i = 0; i < num_intersect; i++) graph_file_close(&intersect_gfiles[i]);

    // Reset graph info
    for(i = 0; i < db_graph.num_of_cols; i++)
      graph_info_init(&db_graph.ginfo[i]);

    // Zero covgs
    memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));

    // Resize graph
    db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, use_ncols, use_ncols, kmers_in_hash);
    intersect_edges = db_graph.col_edges;
    db_graph.col_edges += db_graph.ht.capacity;
  }

  bool kmers_loaded = take_intersect, colours_loaded = false;

  graph_files_merge_mkhdr(out_ctx_path, gfiles, num_gfiles,
                          kmers_loaded, colours_loaded, intersect_edges,
                          intsct_gname_ptr, &db_graph);

  for(i = 0; i < num_gfiles; i++) graph_file_close(&gfiles[i]);

  strbuf_dealloc(&intersect_gname);
  ctx_free(intersect_paths);
  ctx_free(intersect_gfiles);
  ctx_free(gfiles);

  ctx_free(db_graph.col_edges);
  ctx_free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
