#include "global.h"
#include "commands.h"
#include "file_util.h"
#include "util.h"
#include "db_graph.h"
#include "graph_file_reader.h"
#include "graph_format.h"
#include "breakpoint_caller.h"
#include "kmer_occur.h"
#include "seq_reader.h"
#include "gpath_reader.h"
#include "gpath_checks.h"

const char breakpoints_usage[] =
"usage: "CMD" breakpoints [options] <in.ctx> [in2.ctx ..]\n"
"\n"
"  Use trusted assembled genome to call large events.  Output is gzipped.\n"
"  Memory (bytes) is roughly: num_kmers*(8+8+1) + num_ref_kmers*8 + links_mem\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>       Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>    Load path file (can specify multiple times)\n"
"  -o, --out <out.txt.gz>  Save calls (gzipped output) [default: STDOUT]\n"
"  -s, --seq <in>          Trusted input (can specify multiple times)\n"
"  -r, --minref <N>        Require <N> kmers at ref breakpoint [default: "QUOTE_VALUE(DEFAULT_MIN_REF_NKMERS)"]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"paths",        required_argument, NULL, 'p'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"seq",          required_argument, NULL, '1'},
  {"seq",          required_argument, NULL, 's'},
  {"minref",       required_argument, NULL, 'r'},
  {NULL, 0, NULL, 0}
};

int ctx_breakpoints(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *output_file = NULL;
  size_t min_ref_flank = 0;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

  seq_file_t *tmp_sfile;
  SeqFilePtrBuffer sfilebuf;
  seq_file_ptr_buf_alloc(&sfilebuf, 16);

  // tmp args
  size_t i;

  // Arg parsing
  char cmd[100], shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_push(&gpfiles, &tmp_gpfile, 1);
        break;
      case 'r': cmd_check(!min_ref_flank, cmd); min_ref_flank = cmd_uint32_nonzero(cmd, optarg); break;
      case 'o': cmd_check(!output_file, cmd); output_file = optarg; break;
      case '1':
      case 's':
        if((tmp_sfile = seq_open(optarg)) == NULL)
          die("Cannot read --seq file %s", optarg);
        seq_file_ptr_buf_add(&sfilebuf, tmp_sfile);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" breakpoints -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;
  if(min_ref_flank == 0) min_ref_flank = DEFAULT_MIN_REF_NKMERS;

  if(sfilebuf.len == 0) cmd_print_usage("Require at least one --seq file");
  if(optind == argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  char **graph_paths = argv + optind;
  size_t num_gfiles = argc - optind;
  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, gpfiles.b, gpfiles.len, -1);

  //
  // Get file sizes of sequence files
  //
  // set to -1 if we cannot calc
  int64_t est_num_bases = seq_est_seq_bases(sfilebuf.b, sfilebuf.len);
  if(est_num_bases < 0) {
    warn("Cannot get file sizes, using pipes");
    est_num_bases = memargs.num_kmers;
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;
  size_t graph_mem, path_mem;
  int64_t max_req_kmers = MAX2(est_num_bases, (int64_t)ctx_max_kmers);
  int64_t sum_req_kmers = est_num_bases + ctx_sum_kmers;

  // DEV: use threads in memory calculation

  // kmer memory = Edges + paths + 1 bit per colour for in-colour
  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 +
                  (gpfiles.len > 0 ? sizeof(GPath*)*8 : 0) +
                  ncols +
                  sizeof(KONodeList) + sizeof(KOccur) + // see kmer_occur.h
                  8; // 1 byte per kmer for each base to load sequence files

  // For k=31, 8+1+12+8+1=30 bytes per kmer, could be 8+1+8+8=26
  // for human ~3 billion kmers, hash table load factor of 0.75
  //   3*10^9*30 / 0.75 bytes in GiB = 112GB
  //   3*10^9*26 / 0.75 bytes in GiB = 97GB
  // or even (5 instead of 8 bytes for count):
  //   3*10^9*23 / 0.75 bytes in GiB = 86GB

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        max_req_kmers, sum_req_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.b, gpfiles.len, ncols, rem_mem, false);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Open output file
  //
  gzFile gzout = futil_gzopen_create(output_file != NULL ? output_file : "-", "w");

  //
  // Set up memory
  //
  size_t kmer_size = gfiles[0].hdr.kmer_size;

  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL | DBG_ALLOC_BKTLOCKS);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.b, gpfiles.len,
                             path_mem, false, &db_graph);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_load(&gpfiles.b[i], true, &db_graph);

  // Get array of sequence file paths
  size_t num_seq_paths = sfilebuf.len;
  char **seq_paths = ctx_calloc(num_seq_paths, sizeof(char*));
  for(i = 0; i < num_seq_paths; i++)
    seq_paths[i] = strdup(sfilebuf.b[i]->path);

  //
  // Load reference sequence into a read buffer
  //
  ReadBuffer rbuf;
  read_buf_alloc(&rbuf, 1024);
  seq_load_all_reads(sfilebuf.b, sfilebuf.len, &rbuf);

  // Remove commas and colons from read names so we can print:
  //   chr1:start1-end1,chr2:start2-end2...
  for(i = 0; i < rbuf.len; i++) {
    read_t *r = &rbuf.b[i];
    seq_read_truncate_name(r); // strip fast[aq] comments (after whitespace)
    string_char_replace(r->name.b, ',', '.'); // change , -> . in read name
    string_char_replace(r->name.b, ':', ';'); // change : -> ; in read name
  }

  // Create array of cJSON** from input files
  cJSON **hdrs = ctx_malloc(gpfiles.len * sizeof(cJSON*));
  for(i = 0; i < gpfiles.len; i++) hdrs[i] = gpfiles.b[i].json;

  // Call breakpoints
  breakpoints_call(nthreads,
                   gzout, output_file,
                   rbuf.b, rbuf.len,
                   seq_paths, num_seq_paths,
                   min_ref_flank,
                   hdrs, gpfiles.len,
                   &db_graph);

  // Finished: do clean up
  gzclose(gzout);
  ctx_free(hdrs);

  // Close input files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_close(&gpfiles.b[i]);
  gpfile_buf_dealloc(&gpfiles);

  for(i = 0; i < rbuf.len; i++) seq_read_dealloc(&rbuf.b[i]);
  read_buf_dealloc(&rbuf);

  seq_file_ptr_buf_dealloc(&sfilebuf);

  for(i = 0; i < num_seq_paths; i++) free(seq_paths[i]);
  ctx_free(seq_paths);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
