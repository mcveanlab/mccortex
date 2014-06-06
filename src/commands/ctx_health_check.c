#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "graph_file_reader.h"
#include "path_format.h"
#include "path_file_reader.h"
#include "graph_paths.h"
#include "path_store.h"

const char health_usage[] =
"usage: "CMD" check [options] <graph.ctx>\n"
"  Load a graph into memory along with any path files to check they are valid.\n"
"\n"
"  -h, --help             This help message\n"
"  -m, --memory <mem>     Memory to use\n"
"  -n, --nkmers <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>       Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>   Load path file (can specify multiple times)\n"
//
"  -E, --no-edge-check    Don't check kmer edges\n"
"\n";

// Note: although it seems like we should load path files one at a time and
//       check them, that has the down side of not checking merging code.
//       Therefore we load them all at once, which requires more memory.

static struct option longopts[] =
{
// General options
  {"help",          no_argument,       NULL, 'h'},
  {"memory",        required_argument, NULL, 'm'},
  {"nkmers",        required_argument, NULL, 'n'},
  {"threads",       required_argument, NULL, 't'},
  {"paths",         required_argument, NULL, 'p'},
// command specific
  {"no-edge-check", no_argument,       NULL, 'H'},
  {NULL, 0, NULL, 0}
};

int ctx_health_check(int argc, char **argv)
{
  size_t num_of_threads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  bool do_edge_check = true;

  PathFileReader tmp_pfile;
  PathFileBuffer pfilesbuf;
  pfile_buf_alloc(&pfilesbuf, 8);

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 't':
        if(num_of_threads) die("%s set twice", cmd);
        num_of_threads = cmd_parse_arg_uint32_nonzero(cmd, optarg);
        break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'p':
        tmp_pfile = INIT_PATH_READER;
        path_file_open(&tmp_pfile, optarg, true);
        pfile_buf_add(&pfilesbuf, tmp_pfile);
        break;
      case 'E': if(!do_edge_check) die("%s set twice", cmd); do_edge_check=false; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" supernodes -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(optind+1 != argc)
    cmd_print_usage("Too %s arguments", optind == argc ? "few" : "many");

  char *ctx_path = argv[optind];

  if(!do_edge_check && pfilesbuf.len == 0) {
    cmd_print_usage("-E, --no-edge-check and no path files (-p in.ctp). "
                    "Nothing to check.");
  }

  //
  // Open Graph file
  //
  GraphFileReader gfile = INIT_GRAPH_READER;
  graph_file_open(&gfile, ctx_path, true); // true => errors are fatal
  size_t ncols = graph_file_outncols(&gfile);

  //
  // Open path files
  //
  size_t i, path_max_mem = 0, path_max_usedcols = 0;

  for(i = 0; i < pfilesbuf.len; i++) {
    PathFileReader *pfile = &pfilesbuf.data[i];
    path_max_mem = MAX2(path_max_mem, pfile->hdr.num_path_bytes);
    path_max_usedcols = MAX2(path_max_usedcols, path_file_usedcols(pfile));
  }

  // Check for compatibility between graph files and path files
  graphs_paths_compatible(&gfile, 1, pfilesbuf.data, pfilesbuf.len);

  // Decide on memory
  size_t extra_bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  extra_bits_per_kmer = sizeof(Edges) * ncols * 8 + 1; // edges + in_colour
  kmers_in_hash = cmd_get_kmers_in_hash2(memargs.mem_to_use,
                                         memargs.mem_to_use_set,
                                         memargs.num_kmers,
                                         memargs.num_kmers_set,
                                         extra_bits_per_kmer,
                                         gfile.num_of_kmers, gfile.num_of_kmers,
                                         false, &graph_mem);

  // Paths memory
  path_mem = path_files_mem_required(pfilesbuf.data, pfilesbuf.len, false, false,
                                     path_max_usedcols, 0);
  cmd_print_mem(path_mem, "paths");

  total_mem = path_mem + graph_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, ncols, ncols, kmers_in_hash);

  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity * ncols, sizeof(Edges));
  db_graph.node_in_cols = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity)*ncols, 1);

  // Paths
  if(pfilesbuf.len > 0) {
    path_store_alloc(&db_graph.pstore, path_mem, false,
                     db_graph.ht.capacity, path_max_usedcols);
  }

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  graph_load(&gfile, gprefs, NULL);

  // Load path files (if there are any)
  paths_format_merge(pfilesbuf.data, pfilesbuf.len, false, true,
                     num_of_threads, &db_graph);

  // Close files
  for(i = 0; i < pfilesbuf.len; i++) path_file_close(&pfilesbuf.data[i]);
  pfile_buf_dealloc(&pfilesbuf);

  graph_file_close(&gfile);

  if(do_edge_check)
    db_graph_healthcheck(&db_graph);

  if(pfilesbuf.len)
  {
    GraphPathPairing gp;
    gp_alloc(&gp, ncols);
    for(i = 0; i < ncols; i++) gp.ctxcols[i] = gp.ctpcols[i] = i;

    status("Running path check...");
    // Check data store
    path_store_integrity_check(&db_graph.pstore);

    status("  Tracing reads through the graph...");
    graph_paths_check_all_paths(&gp, &db_graph);

    gp_dealloc(&gp);
  }

  status("All looks good!");

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
