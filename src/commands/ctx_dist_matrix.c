#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graphs_load.h"
#include "gpath_checks.h"

const char dist_matrix_usage[] =
"usage: "CMD" dist [options] <in.ctx> [in2.ctx ...]\n"
"\n"
"  Generate a distance matrix counting kmers shared between colours.\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -f, --force           Overwrite output files\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -o, --out <out.csv>   Ouput matrix, tab separated [defaults to STDOUT]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"force",        no_argument,       NULL, 'f'},
  {"out",          required_argument, NULL, 'o'},
  {NULL, 0, NULL, 0}
};

typedef struct {
  const dBGraph *db_graph;
  uint64_t **matrices;
} DistMatrixThreads;

static bool dist_matrix_thread(hkey_t hkey, size_t threadid, void *arg)
{
  DistMatrixThreads *workers = (DistMatrixThreads*)arg;
  const dBGraph *db_graph = workers->db_graph;
  const size_t ncols = db_graph->num_of_cols;
  uint64_t *matrix = workers->matrices[threadid];
  size_t i, j;
  for(i = 0; i < ncols; i++) {
    if(db_node_has_col(db_graph, hkey, i)) {
      matrix[ncols*i+i]++; // colour with self
      for(j = i+1; j < ncols; j++) {
        if(db_node_has_col(db_graph, hkey, j)) {
          matrix[ncols*i+j]++;
        }
      }
    }
  }

  return false; // keep iterating
}

int ctx_dist_matrix(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;

  char *out_path = NULL;

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
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" dist -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults
  if(!nthreads) nthreads = DEFAULT_NTHREADS;

  if(optind >= argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  const size_t num_gfiles = argc - optind;
  char **graph_paths = argv + optind;
  ctx_assert(num_gfiles > 0);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t i, j, ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, NULL, 0, -1);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +

  bits_per_kmer = sizeof(BinaryKmer)*8 + ncols; // kmer + in colour

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  size_t total_mem = graph_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);


  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 0, kmers_in_hash,
                 DBG_ALLOC_NODE_IN_COL);

  // Allocate thread memory
  uint64_t **matrices = ctx_calloc(nthreads, sizeof(DistMatrixThreads));
  DistMatrixThreads workers = {.db_graph = &db_graph, .matrices = matrices};
  for(i = 0; i < nthreads; i++)
    workers.matrices[i] = ctx_calloc(ncols*ncols, sizeof(uint64_t));

  // Open output file
  // Print to stdout unless --out <out> is specified
  FILE *fout = futil_fopen_create(!out_path ? "-" : out_path, "w");


  //
  // Load graphs
  //
  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);
  gprefs.empty_colours = true;

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, NULL);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Generate matrix
  status("[dist_matrix] Generating matrix between %zu colours with %zu thread%s",
         ncols, nthreads, util_plural_str(nthreads));
  hash_table_iterate(&db_graph.ht, nthreads, dist_matrix_thread, &workers);

  // Merge matrices
  for(i = 1; i < nthreads; i++)
    for(j = 0; j < ncols*ncols; j++)
      matrices[0][j] += matrices[i][j];

  size_t row, col;
  uint64_t *mat = matrices[0];

  // Print matrix
  fprintf(fout, ".");// top left column empty
  for(row = 0; row < ncols; row++) fprintf(fout, "\tcol%zu", row);
  fprintf(fout, "\n");
  for(row = 0; row < ncols; row++) {
    fprintf(fout, "col%zu", row);
    for(col = 0; col < ncols; col++)
      fprintf(fout, "\t%zu", (size_t)mat[ncols*row+col]);
    fprintf(fout, "\n");
  }

  status("[dist_matrix]   written to %s", futil_outpath_str(out_path));
  fclose(fout);

  for(i = 0; i < nthreads; i++)
    ctx_free(workers.matrices[i]);
  ctx_free(workers.matrices);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
