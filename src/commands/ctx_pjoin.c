#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
#include "gpath_save.h"

const char pjoin_usage[] =
"usage: "CMD" pjoin [options] <in1.ctp.gz> [[offset:]in2.ctp[:0,2-4] ...]\n"
"\n"
"  Merge cortex path files.\n"
"\n"
"  -h, --help             This help message\n"
"  -q, --quiet            Silence status output normally printed to STDERR\n"
"  -f, --force            Overwrite output files\n"
"  -o, --out <out.ctp.gz> Output file [required]\n"
"  -m, --memory <mem>     Memory to use (required) recommend 80G for human\n"
"  -n, --nkmers <nkmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>      Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
//
"  -g, --graph <in.ctx>   Get number of hash table entries from graph file\n"
"  -c, --outcols <C>      How many 'colours' should the output file have\n"
"  -r, --noredundant      Remove redundant paths\n"
"\n"
"  Files can be specified with specific colours: samples.ctp:2,3\n"
"  Offset specifies where to load the first colour: 3:samples.ctp\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
// command specific
  {"graph",        required_argument, NULL, 'g'},
  {"outcols",      required_argument, NULL, 'c'},
  {"noredundant",  required_argument, NULL, 'r'},
  {NULL, 0, NULL, 0}
};

int ctx_pjoin(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  bool noredundant = false;
  size_t output_ncols = 0;
  char *graph_file = NULL;
  const char *out_ctp_path = NULL;

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
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!out_ctp_path, cmd); out_ctp_path = optarg; break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'g': cmd_check(!graph_file,cmd); graph_file = optarg; break;
      case 'c': cmd_check(!output_ncols, cmd); output_ncols = cmd_uint32_nonzero(cmd, optarg); break;
      case 'r': cmd_check(!noredundant,cmd); noredundant = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" pjoin -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;

  if(out_ctp_path == NULL) cmd_print_usage("--out <out.ctp.gz> required");
  if(optind >= argc) cmd_print_usage("Please specify at least one input file");

  // argi .. argend-1 are graphs to load
  size_t num_pfiles = (size_t)(argc - optind);
  char **paths = argv + optind;

  //
  // Open all path files
  //
  size_t i, j;
  size_t ctp_max_cols = 0;
  uint64_t ctp_max_kmers = 0, ctp_sum_kmers = 0;
  GPathReader *pfiles = ctx_calloc(num_pfiles, sizeof(GPathReader));

  for(i = 0; i < num_pfiles; i++)
  {
    gpath_reader_open2(&pfiles[i], paths[i], "r", ctp_max_cols);

    size_t nkmers = gpath_reader_get_num_kmers(&pfiles[i]);

    ctp_max_cols = MAX2(ctp_max_cols, file_filter_into_ncols(&pfiles[i].fltr));
    ctp_max_kmers = MAX2(ctp_max_kmers, nkmers);
    ctp_sum_kmers += nkmers;

    file_filter_status(&pfiles[i].fltr);
  }


  if(output_ncols == 0) output_ncols = ctp_max_cols;
  else if(ctp_max_cols > output_ncols) {
    cmd_print_usage("You specified --outcols %zu but inputs need at %zu colours",
                    output_ncols, ctp_max_cols);
  }

  // Open graph file to get number of kmers is passed
  GraphFileReader gfile;
  memset(&gfile, 0, sizeof(GraphFileReader));

  if(memargs.num_kmers_set) {
    ctp_max_kmers = MAX2(ctp_max_kmers, memargs.num_kmers);
  }

  if(graph_file != NULL) {
    graph_file_open(&gfile, graph_file);
    ctp_sum_kmers = MIN2(ctp_sum_kmers, graph_file_nkmers(&gfile));
  }

  // Check for compatibility between graph files and path files
  graphs_gpaths_compatible(&gfile, graph_file ? 1 : 0, pfiles, num_pfiles, -1);

  // Done with the graph file now
  if(graph_file != NULL)
    graph_file_close(&gfile);

  if(memargs.num_kmers_set && memargs.num_kmers > ctp_sum_kmers) {
    char num_kmers_str[100], args_num_kmers_str[100];
    ulong_to_str(ctp_sum_kmers, num_kmers_str);
    ulong_to_str(memargs.num_kmers, args_num_kmers_str);
    warn("Using %s kmers instead of (-n) %s", num_kmers_str, args_num_kmers_str);
  }

  // if(num_kmers < ctp_max_kmers) {
  //   cmd_print_usage("Please set a larger -n <kmers> (needs to be > %zu)",
  //               ctp_max_kmers);
  // }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  // Each kmer stores a pointer to its list of paths
  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(GPath*)*8;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctp_max_kmers, ctp_sum_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(pfiles, num_pfiles, output_ncols, rem_mem, true);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  total_mem = graph_mem + path_mem;

  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  // Open output file
  gzFile gzout = futil_gzopen_create(out_ctp_path, "w");

  // Set up graph and PathStore
  size_t kmer_size = gpath_reader_get_kmer_size(&pfiles[0]);
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, output_ncols, 0, kmers_in_hash, 0);

  // Create a path store that tracks path counts
  gpath_reader_alloc_gpstore(pfiles, num_pfiles,
                             path_mem, true, &db_graph);

  for(i = 0; i < num_pfiles; i++)
    gpath_reader_load_sample_names(&pfiles[i], &db_graph);

  // Load contig hist distribution
  ZeroSizeBuffer *contig_histgrms = ctx_calloc(output_ncols, sizeof(ZeroSizeBuffer));

  for(i = 0; i < output_ncols; i++)
    zsize_buf_alloc(&contig_histgrms[i], 512);

  size_t fromcol, intocol;
  for(i = 0; i < num_pfiles; i++) {
    for(j = 0; j < file_filter_num(&pfiles[i].fltr); j++) {
      fromcol = file_filter_fromcol(&pfiles[i].fltr, j);
      intocol = file_filter_intocol(&pfiles[i].fltr, j);
      gpath_reader_load_contig_hist(pfiles[i].json, pfiles[i].fltr.path.b,
                                    fromcol, &contig_histgrms[intocol]);
    }
  }

  // Load path files
  for(i = 0; i < num_pfiles; i++)
    gpath_reader_load(&pfiles[i], GPATH_ADD_MISSING_KMERS, &db_graph);

  status("Got %zu path bytes", (size_t)db_graph.gpstore.path_bytes);

  size_t output_threads = MIN2(nthreads, MAX_IO_THREADS);

  cJSON **hdrs = ctx_calloc(num_pfiles, sizeof(cJSON*));
  for(i = 0; i < num_pfiles; i++) hdrs[i] = pfiles[i].json;

  // Write output file
  gpath_save(gzout, out_ctp_path, output_threads, false,
             NULL, NULL, hdrs, num_pfiles,
             contig_histgrms, output_ncols,
             &db_graph);

  for(i = 0; i < output_ncols; i++)
    zsize_buf_dealloc(&contig_histgrms[i]);

  ctx_free(contig_histgrms);

  gzclose(gzout);
  ctx_free(hdrs);

  // Close ctp files
  // Don't close until now since we were using their headers in the output file
  for(i = 0; i < num_pfiles; i++) gpath_reader_close(&pfiles[i]);
  ctx_free(pfiles);

  char pnum_str[100], pbytes_str[100], pkmers_str[100];
  ulong_to_str(db_graph.gpstore.num_paths, pnum_str);
  bytes_to_str(db_graph.gpstore.path_bytes, 1, pbytes_str);
  ulong_to_str(db_graph.gpstore.num_kmers_with_paths, pkmers_str);

  status("Paths written to: %s\n", out_ctp_path);
  status("  %s paths, %s path-bytes, %s kmers", pnum_str, pbytes_str, pkmers_str);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
