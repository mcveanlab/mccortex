#include "global.h"
#include "commands.h"
#include "range.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
#include "bubble_caller.h"

// Long flanks help us map calls
// increasing allele length can be costly
#define DEFAULT_MAX_FLANK 500
#define DEFAULT_MAX_ALLELE 2000

#define SUBCMD "bubbles"

const char bubbles_usage[] =
"usage: "CMD" "SUBCMD" [options] <in.ctx> [in2.ctx ...]\n"
"\n"
"  Find bubbles in the graph, which are potential variants.\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -o, --out <bub.txt.gz>  Output file [required]\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>       Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>    Load path file (can specify multiple times)\n"
//
"  -H, --haploid <col>     List of haploid colours (e.g. ref colour); '*' means all\n"
"  -A, --max-allele <len>  Max bubble branch length in kmers [default: "QUOTE_VALUE(DEFAULT_MAX_ALLELE)"]\n"
"  -F, --max-flank <len>   Max flank length in kmers [default: "QUOTE_VALUE(DEFAULT_MAX_FLANK)"]\n"
"  -S, --keep-serial       Keep serial bubbles. Use if mapping is hard. Higher FP.\n"
"\n"
"  When loading path files with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into.\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"paths",        required_argument, NULL, 'p'},
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"haploid",      required_argument, NULL, 'H'},
  {"max-allele",   required_argument, NULL, 'A'},
  {"max-flank",    required_argument, NULL, 'F'},
  {"keep-serial",  required_argument, NULL, 'S'},
  {NULL, 0, NULL, 0}
};

int ctx_bubbles(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t max_allele_len = 0, max_flank_len = 0;
  bool remove_serial_bubbles = true;

  // List of haploid colours
  size_t *hapcols = NULL;
  int nhapcols = 0;
  char *hapcols_arg = NULL;

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
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_push(&gpfiles, &tmp_gpfile, 1);
        break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'H': cmd_check(!hapcols_arg, cmd); hapcols_arg = optarg; break;
      case 'A': cmd_check(!max_allele_len, cmd); max_allele_len = cmd_uint32_nonzero(cmd, optarg); break;
      case 'F': cmd_check(!max_flank_len, cmd); max_flank_len = cmd_uint32_nonzero(cmd, optarg); break;
      case 'S': cmd_check(remove_serial_bubbles,cmd); remove_serial_bubbles = false; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" "SUBCMD" -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;
  if(max_allele_len == 0) max_allele_len = DEFAULT_MAX_ALLELE;
  if(max_flank_len == 0) max_flank_len = DEFAULT_MAX_FLANK;

  if(optind >= argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  const size_t num_gfiles = argc - optind;
  char **graph_paths = argv + optind;
  ctx_assert(num_gfiles > 0);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t i, ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, gpfiles.b, gpfiles.len, -1);

  //
  // Check haploid colours are valid
  //
  if(hapcols_arg != NULL) {
    if((nhapcols = range_get_num(hapcols_arg, ncols)) < 0)
      die("Invalid haploid colour list: %s", hapcols_arg);

    hapcols = ctx_calloc(nhapcols, sizeof(hapcols[0]));
    if(range_parse_array(hapcols_arg, hapcols, ncols) < 0)
      die("Invalid haploid colour list: %s", hapcols_arg);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, thread_mem;
  char thread_mem_str[100];

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +
  // visitedfw/rv(2bits/thread)

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 +
                  (gpfiles.len > 0 ? sizeof(GPath*)*8 : 0) +
                  ncols + 2*nthreads;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, &graph_mem);

  // Thread memory
  thread_mem = roundup_bits2bytes(kmers_in_hash) * 2;
  bytes_to_str(thread_mem * nthreads, 1, thread_mem_str);
  status("[memory] (of which threads: %zu x %zu = %s)\n",
          nthreads, thread_mem, thread_mem_str);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem+thread_mem);
  path_mem = gpath_reader_mem_req(gpfiles.b, gpfiles.len, ncols, rem_mem, false);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  size_t total_mem = graph_mem + thread_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Open output file
  //
  gzFile gzout = futil_gzopen_create(out_path, "w");

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.b, gpfiles.len, path_mem, false, &db_graph);

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

  // Load path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_load(&gpfiles.b[i], GPATH_DIE_MISSING_KMERS, &db_graph);

  // Create array of cJSON** from input files
  cJSON **hdrs = ctx_malloc(gpfiles.len * sizeof(cJSON*));
  for(i = 0; i < gpfiles.len; i++) hdrs[i] = gpfiles.b[i].json;

  // Now call variants
  BubbleCallingPrefs call_prefs = {.max_allele_len = max_allele_len,
                                   .max_flank_len = max_flank_len,
                                   .haploid_cols = hapcols,
                                   .nhaploid_cols = nhapcols,
                                   .remove_serial_bubbles = remove_serial_bubbles};

  invoke_bubble_caller(nthreads, &call_prefs,
                       gzout, out_path,
                       hdrs, gpfiles.len,
                       &db_graph);

  status("  saved to: %s\n", out_path);
  gzclose(gzout);
  ctx_free(hdrs);

  // Close input path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_close(&gpfiles.b[i]);
  gpfile_buf_dealloc(&gpfiles);

  ctx_free(hapcols);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
