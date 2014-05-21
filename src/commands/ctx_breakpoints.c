#include "global.h"
#include "commands.h"
#include "file_util.h"
#include "util.h"
#include "db_graph.h"
#include "graph_file_filter.h"
#include "path_file_filter.h"
#include "graph_paths.h"
#include "path_format.h"
#include "path_store.h"
#include "breakpoint_caller.h"
#include "kmer_occur.h"
#include "seq_reader.h"

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
"  Use trusted assembled genome to call large events.  Output is gzipped.\n"
"\n"
"  --seq <in>          Trusted input (can be used multiple times)\n"
"  --out <out.txt.gz>  Save calls (gzipped output) [default: STDOUT]\n"
"  --minref <N>        Require <N> kmers at ref breakpoint [default: "QUOTE_VALUE(DEFAULT_MIN_REF_NKMERS)"]\n"
"  --maxref <N>        Limit to <N> kmers at ref breakpoint [default: "QUOTE_VALUE(DEFAULT_MAX_REF_NKMERS)"]\n";

int ctx_breakpoints(CmdArgs *args)
{
  int argi, argc = args->argc;
  char **argv = args->argv;
  // Already checked that input has at least 1 argument

  size_t max_seq_inputs = argc / 2;
  seq_file_t **seq_files = ctx_malloc(max_seq_inputs * sizeof(seq_file_t*));
  const char **seq_paths = ctx_malloc(max_seq_inputs * sizeof(char*));
  size_t i, num_seq_files = 0;
  size_t min_ref_flank = DEFAULT_MIN_REF_NKMERS;
  size_t max_ref_flank = DEFAULT_MAX_REF_NKMERS;
  size_t est_num_bases = 0;

  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++) {
    if(strcmp(argv[argi],"--seq") == 0) {
      if(argi+1 == argc) cmd_print_usage("--seq <in.fa> requires a file");

      char *file_path = argv[argi+1];
      if((seq_files[num_seq_files] = seq_open(file_path)) == NULL)
        die("Cannot read --seq file %s", file_path);

      off_t fsize = futil_get_file_size(file_path);
      if(fsize < 0) warn("Cannot get file size: %s", file_path);
      else {
        if(seq_is_fastq(seq_files[num_seq_files]) ||
           seq_is_sam(seq_files[num_seq_files])) {
          est_num_bases += fsize / 2;
        } else {
          est_num_bases += fsize;
        }
      }

      seq_paths[num_seq_files] = file_path;
      num_seq_files++;
      argi++; // We took an argument
    }
    else if(!strcmp(argv[argi],"--minref")) {
      if(argi+1 == argc) die("%s <N> requires an argument", argv[argi]);
      if(!parse_entire_size(argv[argi+1], &min_ref_flank) || min_ref_flank==0) {
        die("Invalid argument %s %s must be >= 1", argv[argi], argv[argi+1]);
      }
      argi++; // We took an argument
    }
    else if(!strcmp(argv[argi],"--maxref")) {
      if(argi+1 == argc) die("%s <N> requires an argument", argv[argi]);
      if(!parse_entire_size(argv[argi+1], &max_ref_flank) || max_ref_flank==0) {
        die("Invalid argument %s %s must be >= 1", argv[argi], argv[argi+1]);
      }
      argi++; // We took an argument
    }
    else {
      cmd_print_usage("Unknown arg: %s", argv[argi]);
    }
  }

  if(min_ref_flank > max_ref_flank) {
    warn("Setting max ref flank len to be min (was: %zu now: %zu)",
         max_ref_flank, min_ref_flank);
    max_ref_flank = min_ref_flank;
  }

  if(num_seq_files == 0) cmd_print_usage("Require at least one --seq file");
  if(argi == argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  char **graph_paths = argv + argi;
  size_t num_gfiles = argc - argi;
  GraphFileReader gfiles[num_gfiles];
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Open path files
  //
  size_t num_pfiles = args->num_ctp_files;
  PathFileReader pfiles[num_pfiles];
  size_t path_max_usedcols = 0;

  for(i = 0; i < num_pfiles; i++) {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], args->ctp_files[i], true);
    path_max_usedcols = MAX2(path_max_usedcols, path_file_usedcols(&pfiles[i]));
  }

  // Check graph + paths are compatible
  graphs_paths_compatible(gfiles, num_gfiles, pfiles, num_pfiles);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;
  size_t graph_mem, path_mem;
  size_t max_req_kmers = MAX2(est_num_bases, ctx_max_kmers);
  size_t sum_req_kmers = est_num_bases + ctx_sum_kmers;

  // DEV: use threads in memory calculation
  size_t num_of_threads = args->max_work_threads;

  // DEV: pass path memory into cmd_get_kmers_in_hash
  // Path memory
  path_mem = path_files_mem_required(pfiles, num_pfiles, false, false);
  cmd_print_mem(path_mem, "paths");

  // kmer memory = Edges + paths + 1 bit per colour for in-colour
  bits_per_kmer = sizeof(Edges)*8 + sizeof(PathIndex)*8 + ncols +
                  sizeof(KONodeList) + sizeof(KOccur) + // see kmer_occur.h
                  8; // 1 byte per kmer for each base to load sequence files

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        max_req_kmers, sum_req_kmers,
                                        false, &graph_mem);

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  //
  // Open output file
  //
  const char *output_file = args->output_file_set ? args->output_file : "-";
  gzFile gzout = futil_gzopen_output(output_file);

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
  path_store_alloc(&db_graph.pstore, path_mem, false,
                   db_graph.ht.capacity, path_max_usedcols);

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

  hash_table_print_stats(&db_graph.ht);

  //
  // Load path files (does nothing if num_pfiles == 0)
  paths_format_merge(pfiles, num_pfiles, false, false, num_of_threads, &db_graph);

  for(i = 0; i < num_pfiles; i++) path_file_close(&pfiles[i]);

  //
  // Load reference sequence into a read buffer
  //
  ReadBuffer rbuf;
  readbuf_alloc(&rbuf, 1024);
  seq_load_all_reads(seq_files, num_seq_files, &rbuf);

  // Remove commas and colons from read names so we can print:
  //   chr1:start1-end1,chr2:start2-end2...
  for(i = 0; i < rbuf.len; i++) {
    read_t *r = &rbuf.data[i];
    seq_read_truncate_name(r); // strip fast[aq] comments (after whitespace)
    string_char_replace(r->name.b, ',', '.'); // change , -> . in read name
    string_char_replace(r->name.b, ':', ';'); // change : -> ; in read name
  }

  // Call breakpoints
  breakpoints_call(num_of_threads, rbuf.data, rbuf.len,
                   gzout, output_file, seq_paths, num_seq_files,
                   min_ref_flank, max_ref_flank,
                   args, &db_graph);

  // Finished: do clean up
  gzclose(gzout);

  for(i = 0; i < rbuf.len; i++) seq_read_dealloc(&rbuf.data[i]);
  readbuf_dealloc(&rbuf);

  ctx_free(seq_files);
  ctx_free(seq_paths);
  ctx_free(db_graph.col_edges);
  ctx_free(db_graph.node_in_cols);
  ctx_free(db_graph.bktlocks);

  path_store_dealloc(&db_graph.pstore);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
