#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "gpath_reader.h"
#include "gpath_checks.h"

const char health_usage[] =
"usage: "CMD" check [options] <graph.ctx>\n"
"  Load a graph into memory along with any path files to check they are valid.\n"
"\n"
"  -h, --help             This help message\n"
"  -q, --quiet            Silence status output normally printed to STDERR\n"
"  -m, --memory <mem>     Memory to use\n"
"  -n, --nkmers <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>      Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
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
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  bool do_edge_check = true;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

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
        if(nthreads) die("%s set twice", cmd);
        nthreads = cmd_uint32_nonzero(cmd, optarg);
        break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_push(&gpfiles, &tmp_gpfile, 1);
        break;
      case 'E': if(!do_edge_check) die("%s set twice", cmd); do_edge_check=false; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" check -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;

  if(optind+1 != argc)
    cmd_print_usage("Too %s arguments", optind == argc ? "few" : "many");

  const char *ctx_path = argv[optind];

  if(!do_edge_check && gpfiles.len == 0) {
    cmd_print_usage("-E, --no-edge-check and no path files (-p in.ctp). "
                    "Nothing to check.");
  }

  //
  // Open Graph file
  //
  GraphFileReader gfile;
  memset(&gfile, 0, sizeof(GraphFileReader));
  graph_file_open(&gfile, ctx_path);
  size_t ncols = file_filter_into_ncols(&gfile.fltr);

  // Check for compatibility between graph files and path files
  graphs_gpaths_compatible(&gfile, 1, gpfiles.b, gpfiles.len, -1);

  //
  // Decide on memory
  //
  size_t i, bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  // edges + in_colour
  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges) * ncols * 8 + 1 +
                  (gpfiles.len > 0 ? sizeof(GPath*)*8 : 0);

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        gfile.num_of_kmers, gfile.num_of_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.b, gpfiles.len, ncols, rem_mem, false);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  total_mem = path_mem + graph_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, ncols, ncols, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.b, gpfiles.len, path_mem, false, &db_graph);

  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);
  gprefs.empty_colours = true;

  graph_load(&gfile, gprefs, NULL);

  // Load path files
  for(i = 0; i < gpfiles.len; i++) {
    gpath_reader_load(&gpfiles.b[i], GPATH_DIE_MISSING_KMERS, &db_graph);
    gpath_reader_close(&gpfiles.b[i]);
  }

  if(do_edge_check)
    db_graph_healthcheck(&db_graph);

  if(gpfiles.len) {
    status("Tracing reads through the graph...");
    gpath_checks_all_paths(&db_graph, nthreads);
  }

  graph_file_close(&gfile);
  gpfile_buf_dealloc(&gpfiles);
  db_graph_dealloc(&db_graph);

  status("All looks good!");

  return EXIT_SUCCESS;
}
