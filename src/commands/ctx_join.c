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
#include "graph_file_reader.h"

// Given (A,B,C) are ctx binaries, A:1 means colour 1 in A,
// {A:1,B:0} is loading A:1 and B:0 into a single colour
//
// Default behaviour is to load colours consecutively
//   output: {A:0},{A:1},{B:0},{B:1},{B:2},{C:0},{C:1}
//

const char join_usage[] =
"usage: "CMD" join [options] in1.ctx [[offset:]in2.ctx[:1,2,4-5] ...]\n"
"\n"
"  Merge cortex graphs.\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -o, --out <out.ctx>     Output file [required]\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
//
"  -N, --ncols <c>         How many colours to load at once [default: 1]\n"
"  -i, --intersect <a.ctx> Only load the kmers that are in graph A.ctx. Can be\n"
"                          specified multiple times. <a.ctx> is NOT merged into\n"
"                          the output file.\n"
"\n"
"  Files can be specified with specific colours: samples.ctx:2,3\n"
"  Offset specifies where to load the first colour: 3:samples.ctx\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
// command specific
  {"ncols",        required_argument, NULL, 'N'},
  {"intersect",    required_argument, NULL, 'i'},
  {NULL, 0, NULL, 0}
};

static inline void remove_non_intersect_nodes(hkey_t node, Covg *covgs,
                                              Covg num, HashTable *ht)
{
  if(covgs[node] != num)
    hash_table_delete(ht, node);
}

int ctx_join(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t use_ncols = 0;

  GraphFileReader tmp_gfile;
  GraphFileBuffer isec_gfiles_buf;
  gfile_buf_alloc(&isec_gfiles_buf, 8);

  // Arg parsing
  char cmd[100], shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'N': cmd_check(!use_ncols, cmd); use_ncols = cmd_uint32_nonzero(cmd, optarg); break;
      case 'i':
        memset(&tmp_gfile, 0, sizeof(GraphFileReader));
        graph_file_open(&tmp_gfile, optarg);
        if(!file_filter_into_ncols(&tmp_gfile.fltr) > 1)
          warn("Flattening intersection graph into colour 0: %s", optarg);
        file_filter_flatten(&tmp_gfile.fltr, 0);
        gfile_buf_add(&isec_gfiles_buf, tmp_gfile);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" join -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  GraphFileReader *igfiles = isec_gfiles_buf.data;
  size_t num_igfiles = isec_gfiles_buf.len;

  bool use_ncols_set = (use_ncols > 0);

  if(use_ncols == 0) use_ncols = 1;
  if(!out_path) cmd_print_usage("--out <out.ctx> required");

  if(optind >= argc)
    cmd_print_usage("Please specify at least one input graph file");

  // optind .. argend-1 are graphs to load
  size_t num_gfiles = (size_t)(argc - optind);
  char **gfile_paths = argv + optind;

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));

  status("Probing %zu graph files and %zu intersect files", num_gfiles, num_igfiles);

  // Check all binaries are valid binaries with matching kmer size
  size_t i;
  size_t ctx_max_cols = 0;
  uint64_t min_intersect_num_kmers = 0, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    graph_file_open2(&gfiles[i], gfile_paths[i], "r", ctx_max_cols);

    if(gfiles[0].hdr.kmer_size != gfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                      gfiles[0].hdr.kmer_size, gfiles[i].hdr.kmer_size);
    }

    ctx_max_cols = MAX2(ctx_max_cols, file_filter_into_ncols(&gfiles[i].fltr));
    ctx_max_kmers = MAX2(ctx_max_kmers, graph_file_nkmers(&gfiles[i]));
    ctx_sum_kmers += gfiles[i].num_of_kmers;
  }

  // Probe intersection graph files
  for(i = 0; i < num_igfiles; i++)
  {
    if(gfiles[0].hdr.kmer_size != igfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                  gfiles[0].hdr.kmer_size, igfiles[i].hdr.kmer_size);
    }

    uint64_t nkmers = graph_file_nkmers(&igfiles[i]);

    if(i == 0) min_intersect_num_kmers = nkmers;
    else if(nkmers < min_intersect_num_kmers)
    {
      // Put smallest intersection binary first
      SWAP(igfiles[i], igfiles[0]);
      min_intersect_num_kmers = nkmers;
    }
  }

  bool take_intersect = (num_igfiles > 0);

  if(take_intersect)
    ctx_max_kmers = min_intersect_num_kmers;

  if(use_ncols < ctx_max_cols && strcmp(out_path,"-") == 0)
    die("I need %zu colours if outputting to STDOUT (--ncols)", ctx_max_cols);

  // Check out_path is writable
  futil_create_output(out_path);

  if(use_ncols > ctx_max_cols) {
    warn("I only need %zu colour%s ('--ncols %zu' ignored)",
         ctx_max_cols, util_plural_str(ctx_max_cols), use_ncols);
    use_ncols = ctx_max_cols;
  }

  status("Output %zu cols; from %zu files; intersecting %zu graphs; ",
         ctx_max_cols, num_gfiles, num_igfiles);

  if(num_gfiles == 1 && num_igfiles == 0)
  {
    // Loading only one file with no intersection files
    // don't need to store a graph in memory, can filter as stream
    graph_stream_filter_mkhdr(out_path, &gfiles[0], NULL, NULL, NULL);
    graph_file_close(&gfiles[0]);
    gfile_buf_dealloc(&isec_gfiles_buf);
    ctx_free(gfiles);

    return EXIT_SUCCESS;
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  (sizeof(Covg) + sizeof(Edges)) * 8 * use_ncols;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  if(!use_ncols_set)
  {
    // Maximise use_ncols
    size_t max_usencols = (memargs.mem_to_use*8) / bits_per_kmer;
    use_ncols = MIN2(max_usencols, ctx_max_cols);

    // Re-check memory used
    kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                          memargs.mem_to_use_set,
                                          memargs.num_kmers,
                                          memargs.num_kmers_set,
                                          bits_per_kmer,
                                          ctx_max_kmers, ctx_sum_kmers,
                                          true, &graph_mem);
  }

  status("Using %zu colour%s in memory", use_ncols, util_plural_str(use_ncols));

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  // Create db_graph
  dBGraph db_graph;
  Edges *intersect_edges = NULL;
  size_t edge_cols = (use_ncols + take_intersect);

  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, use_ncols, use_ncols,
                 kmers_in_hash, DBG_ALLOC_COVGS);

  // We allocate edges ourself since it's a special case
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity*edge_cols, sizeof(Edges));

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

    for(i = 0; i < num_igfiles; i++)
    {
      graph_load(&igfiles[i], gprefs, NULL);

      // Update intersect header
      // note: intersection graphs all load exactly one colour into colour 0
      graph_info_make_intersect(&igfiles[i].hdr.ginfo[0], &intersect_gname);

      gprefs.must_exist_in_graph = true;
      gprefs.must_exist_in_edges = db_graph.col_edges;
    }

    if(num_igfiles > 1)
    {
      // Remove nodes where covg != num_igfiles
      HASH_ITERATE_SAFE(&db_graph.ht, remove_non_intersect_nodes,
                        db_graph.col_covgs, (Covg)num_igfiles, &db_graph.ht);
    }

    status("Loaded intersection set\n");
    intsct_gname_ptr = intersect_gname.b;

    for(i = 0; i < num_igfiles; i++) graph_file_close(&igfiles[i]);

    // Reset graph info
    for(i = 0; i < db_graph.num_of_cols; i++)
      graph_info_init(&db_graph.ginfo[i]);

    // Zero covgs
    memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));

    // Use union edges we loaded to intersect new edges
    intersect_edges = db_graph.col_edges;
    db_graph.col_edges += db_graph.ht.capacity;
  }

  bool kmers_loaded = take_intersect, colours_loaded = false;

  graph_files_merge_mkhdr(out_path, gfiles, num_gfiles,
                          kmers_loaded, colours_loaded, intersect_edges,
                          intsct_gname_ptr, &db_graph);

  if(take_intersect)
    db_graph.col_edges -= db_graph.ht.capacity;

  for(i = 0; i < num_gfiles; i++) graph_file_close(&gfiles[i]);

  strbuf_dealloc(&intersect_gname);
  gfile_buf_dealloc(&isec_gfiles_buf);
  ctx_free(gfiles);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
