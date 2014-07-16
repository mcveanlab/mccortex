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

/*
 Output format:

##fileFormat=CtxBreakpointsv0.1
##fileDate=20140515
##cmd="../..//bin/ctx31 breakpoints -t 1 -m 10M --seq seq0.fa --seq seq1.fa --out breakpoints.txt.gz mix.k11.ctx"
##wkdir=/Users/isaac/ninja-cortex/test/breakpoint
##reference=seq0.fa:seq1.fa
##ctxVersion="ctx=b871a64 zlib=1.2.5 htslib=0.2.0-rc8-6-gd49dfa6 ASSERTS=ON CHECKS=ON k=3..31"
##ctxKmerSize=31
>brkpnt.0.5pflank chr=seq0:1-20:+:1
CCCGTAGGTAAGGGCGTTAG
>brkpnt.0.3pflank chr=seq1:21-50:+:11
CGGGTTGGAGTTGGCCAAAGAAGTTCAACG
>brkpnt.0.path cols=0
A

>brkpnt.1.5pflank chr=seq0:1-20:+:1
...

*/

const char breakpoints_usage[] =
"usage: "CMD" breakpoints [options] <in.ctx> [in2.ctx ..]\n"
"\n"
"  Use trusted assembled genome to call large events.  Output is gzipped.\n"
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
"  -R, --maxref <N>        Limit to <N> kmers at ref breakpoint [default: "QUOTE_VALUE(DEFAULT_MAX_REF_NKMERS)"]\n"
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
  {"maxref",       required_argument, NULL, 'R'},
  {NULL, 0, NULL, 0}
};

int ctx_breakpoints(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *output_file = NULL;
  size_t min_ref_flank = DEFAULT_MIN_REF_NKMERS;
  size_t max_ref_flank = DEFAULT_MAX_REF_NKMERS;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

  seq_file_t *tmp_sfile;
  SeqFilePtrBuffer sfilebuf;
  seq_file_ptr_buf_alloc(&sfilebuf, 16);

  // tmp args
  size_t i;
  size_t set_min_flank = 0, set_max_flank = 0;

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
        gpfile_buf_add(&gpfiles, tmp_gpfile);
        break;
      case 'r': min_ref_flank = cmd_uint32_nonzero(cmd, optarg); set_min_flank++; break;
      case 'R': max_ref_flank = cmd_uint32_nonzero(cmd, optarg); set_max_flank++; break;
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

  if(set_min_flank > 1) cmd_print_usage("Set --minref more than once");
  if(set_max_flank > 1) cmd_print_usage("Set --maxref more than once");

  if(min_ref_flank > max_ref_flank) {
    warn("Setting max ref flank len to be min (was: %zu now: %zu)",
         max_ref_flank, min_ref_flank);
    max_ref_flank = min_ref_flank;
  }

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
  graphs_gpaths_compatible(gfiles, num_gfiles, gpfiles.data, gpfiles.len);

  //
  // Get file sizes of sequence files
  //
  off_t fsize;
  size_t est_num_bases = 0;

  for(i = 0; i < sfilebuf.len; i++) {
    tmp_sfile = sfilebuf.data[i];
    fsize = futil_get_file_size(tmp_sfile->path);
    if(fsize < 0) warn("Cannot get file size: %s", tmp_sfile->path);
    else {
      if(seq_is_fastq(tmp_sfile) || seq_is_sam(tmp_sfile))
        est_num_bases += fsize / 2;
      else
        est_num_bases += fsize;
    }
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;
  size_t graph_mem, path_mem;
  int64_t max_req_kmers = MAX2(est_num_bases, ctx_max_kmers);
  int64_t sum_req_kmers = est_num_bases + ctx_sum_kmers;

  // DEV: use threads in memory calculation

  // kmer memory = Edges + paths + 1 bit per colour for in-colour
  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 +
                  (gpfiles.len > 0 ? sizeof(GPath*)*8 : 0) +
                  ncols +
                  sizeof(KONodeList) + sizeof(KOccur) + // see kmer_occur.h
                  8; // 1 byte per kmer for each base to load sequence files

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        max_req_kmers, sum_req_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.data, gpfiles.len, ncols, rem_mem, false);

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
  db_graph_alloc(&db_graph, kmer_size, ncols, 1, kmers_in_hash);
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));

  // In colour
  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);
  db_graph.node_in_cols = ctx_calloc(bytes_per_col*ncols, 1);

  // Multithreading locks for the hash table
  db_graph.bktlocks = ctx_calloc(roundup_bits2bytes(db_graph.ht.num_of_buckets), 1);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.data, gpfiles.len,
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
  for(i = 0; i < gpfiles.len; i++) {
    gpath_reader_load(&gpfiles.data[i], true, &db_graph);
    gpath_reader_close(&gpfiles.data[i]);
  }
  gpfile_buf_dealloc(&gpfiles);

  // Get array of sequence file paths
  size_t num_seq_paths = sfilebuf.len;
  char **seq_paths = ctx_calloc(num_seq_paths, sizeof(char*));
  for(i = 0; i < num_seq_paths; i++)
    seq_paths[i] = strdup(sfilebuf.data[i]->path);

  //
  // Load reference sequence into a read buffer
  //
  ReadBuffer rbuf;
  readbuf_alloc(&rbuf, 1024);
  seq_load_all_reads(sfilebuf.data, sfilebuf.len, &rbuf);

  // Remove commas and colons from read names so we can print:
  //   chr1:start1-end1,chr2:start2-end2...
  for(i = 0; i < rbuf.len; i++) {
    read_t *r = &rbuf.data[i];
    seq_read_truncate_name(r); // strip fast[aq] comments (after whitespace)
    string_char_replace(r->name.b, ',', '.'); // change , -> . in read name
    string_char_replace(r->name.b, ':', ';'); // change : -> ; in read name
  }

  // Call breakpoints
  breakpoints_call(nthreads,
                   gzout, output_file,
                   rbuf.data, rbuf.len,
                   seq_paths, num_seq_paths,
                   min_ref_flank, max_ref_flank,
                   &db_graph);

  // Finished: do clean up
  gzclose(gzout);

  for(i = 0; i < rbuf.len; i++) seq_read_dealloc(&rbuf.data[i]);
  readbuf_dealloc(&rbuf);

  seq_file_ptr_buf_dealloc(&sfilebuf);

  for(i = 0; i < num_seq_paths; i++) free(seq_paths[i]);
  ctx_free(seq_paths);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
