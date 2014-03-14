#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "graph_file_filter.h"
#include "path_format.h"
#include "path_file_filter.h"
#include "graph_paths.h"
#include "path_store.h"

const char health_usage[] =
"usage: "CMD" check [options] <graph.ctx>\n"
"  Load a graph into memory along with any path files to check they are valid.\n"
"\n"
"  Options:\n"
"   -m <mem>       Memory to use\n"
"   -n <kmers>     Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   -p <in.ctp>    Load paths\n"
"   --noedgecheck  Don't check kmer edges\n";

// DEV: should load path files one at a time and check them?
//      doing that won't check merging code.

int ctx_health_check(CmdArgs *args)
{
  int argi, argc = args->argc;
  char **argv = args->argv;
  // Have already check that we have exactly 1 argument

  size_t i;
  bool do_edge_check = true;

  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++) {
    if(strcmp(argv[argi],"--noedgecheck") == 0) do_edge_check = false;
    else cmd_print_usage("Unknown option: %s", argv[argi]);
  }

  if(argi+1 < argc) cmd_print_usage("Too many arguments");
  if(argi+1 > argc) cmd_print_usage("Too few arguments");

  char *ctx_path = argv[argi];

  if(!do_edge_check && args->num_ctp_files == 0) {
    cmd_print_usage("--noedgecheck and no path files (-p in.ctp) - nothing to check.");
  }

  //
  // Open Graph file
  //
  GraphFileReader gfile = INIT_GRAPH_READER;
  graph_file_open(&gfile, ctx_path, true); // true => errors are fatal
  size_t col, ncols = graph_file_outncols(&gfile);

  //
  // Open path files
  //
  size_t num_pfiles = args->num_ctp_files;
  PathFileReader pfiles[num_pfiles];
  size_t path_max_mem = 0, path_max_usedcols = 0;

  for(i = 0; i < num_pfiles; i++) {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], args->ctp_files[i], true);
    path_max_mem = MAX2(path_max_mem, pfiles[i].hdr.num_path_bytes);
    path_max_usedcols = MAX2(path_max_usedcols, path_file_usedcols(&pfiles[i]));
  }

  // Check for compatibility between graph files and path files
  graphs_paths_compatible(&gfile, 1, pfiles, num_pfiles);

  // Decide on memory
  size_t extra_bits_per_kmer, kmers_in_hash, graph_mem;
  size_t path_mem = 0, path_mem_req = 0, tmp_path_mem = 0, total_mem;

  extra_bits_per_kmer = sizeof(Edges) * ncols * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        gfile.num_of_kmers, false, &graph_mem);

  // Path Memory
  if(num_pfiles) {
    tmp_path_mem = path_files_tmp_mem_required(pfiles, num_pfiles);
    path_mem_req = path_max_mem + tmp_path_mem;
    // If we don't have enough memory, will bork out at cmd_check_mem_limit
    // Only need the required amount of memory if loading a single paths file
    // May need more if having to merge path files
    if(num_pfiles > 1)
      path_mem = MAX2(args->mem_to_use - graph_mem, path_mem_req);
    else
      path_mem = path_mem_req;

    char path_mem_str[100];
    bytes_to_str(path_mem, 1, path_mem_str);
    status("[memory] paths: %s", path_mem_str);
  }

  total_mem = path_mem + graph_mem;
  cmd_check_mem_limit(args, total_mem);

  // Create db_graph
  dBGraph db_graph;
  size_t num_edge_cols = do_edge_check ? ncols : 1;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, ncols, num_edge_cols, kmers_in_hash);

  // Only need one edge per colour if doing edge check
  if(do_edge_check) {
    db_graph.col_edges = calloc2(db_graph.ht.capacity * ncols, sizeof(Edges));
  } else {
    db_graph.node_in_cols = calloc2(roundup_bits2bytes(db_graph.ht.capacity)*ncols, 1);
    db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  }

  // Paths
  if(num_pfiles > 0) {
    path_store_alloc(&db_graph.pstore, path_mem-tmp_path_mem, tmp_path_mem,
                     db_graph.ht.capacity, path_max_usedcols);
  }

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  graph_load(&gfile, gprefs, NULL);

  // Load path files
  if(num_pfiles) {
    paths_format_merge(pfiles, num_pfiles, false, &db_graph);
  }

  if(do_edge_check) {
    status("Running edge check...");
    db_graph_healthcheck(&db_graph);
  }

  if(num_pfiles) {
    GraphPathPairing gp;
    gp_alloc(&gp, ncols);
    for(i = 0; i < ncols; i++) gp.ctxcols[i] = gp.ctpcols[i] = i;

    status("Running path check...");
    // Check data store
    path_store_integrity_check(&db_graph.pstore);

    status("  Tracing reads through the graph...");
    for(col = 0; col < ncols; col++)
      graph_paths_check_all_paths(&gp, &db_graph);

    gp_dealloc(&gp);
  }

  status("All looks good!");

  if(num_pfiles) path_store_dealloc(&db_graph.pstore);
  if(db_graph.node_in_cols) free(db_graph.node_in_cols);
  free(db_graph.col_edges);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
