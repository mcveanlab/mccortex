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

#include "seq_file.h"

/*
 Output format:

##fileFormat=CtxBreakpointsv0.1
##fileDate=20140204
##reference=file://ref/ref.fa
##cmd="../bin/ctx31 breakpoints --ref ref.fa --out breakpoint.fa thing.ctx"
##wkdir=/home/isaac
##ctxVersion=<version=8d5ff69,MAXK=31>
##ctxKmerSize=31
>101.5pflank colours=1,2
GATTTTTGCACATTCGTGAT [chr1:100-119:+,chr13:181-200:-];
>101.3pflank
GATTTTTGCACATTCGTGAT [chr1:100-119:+,chr13:181-200:-];
>101.path
  GATTTTTGCACATTCGTGAT [chr1:100-119:+,chr13:181-200:-];
  GATTTTTGCACA [];
  GATATTCGTGAT [chr1:100-119:+,chr13:181-200:-];
  TTTTGCACTTCTGAT [chr1:100-119:+,chr13:181-200:-];
>102.5pflank colours=0
GATTTTTGCACATTCGTGAT [chrX:100-119:+];
>102.3pflank
GATTTTTGCACATTCGTGAT [chrY:100-119:+];
>102.path

*/

const char breakpoints_usage[] =
"usage: "CMD" breakpoints [options] <in.ctx> [in2.ctx ..]\n"
"\n"
"  --seq <in>       Trusted input (can be used multiple times)\n"
"  --out <out.txt>  Save calls [default: STDOUT]\n";

#include "objbuf_macro.h"
create_objbuf(readbuf, ReadBuffer, read_t);

int ctx_breakpoints(CmdArgs *args)
{
  int argi, argc = args->argc;
  char **argv = args->argv;
  // Already checked that input has at least 1 argument

  size_t max_seq_inputs = argc / 2;
  seq_file_t *seq_files[max_seq_inputs];
  size_t i, num_seq_files = 0;

  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++) {
    if(strcmp(argv[argi],"--seq") == 0) {
      if(argi+1 == argc) cmd_print_usage("--seq <in.fa> requires a file");
      if((seq_files[num_seq_files] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read --seq file %s", argv[argi+1]);
      num_seq_files++;
      argi++; // We took an argument
    }
    else {
      cmd_print_usage("Unknown arg: %s", argv[argi]);
    }
  }

  if(num_seq_files == 0) cmd_print_usage("Require at least one --seq file");
  if(argi == argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  char **graph_paths = argv + argi;
  size_t num_gfiles = argc - argi;
  GraphFileReader gfiles[num_gfiles];
  size_t ncols = 0, ctx_max_kmers = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    gfiles[i] = INIT_GRAPH_READER;
    graph_file_open(&gfiles[i], graph_paths[i], true);

    // Pile colours on top of each other
    file_filter_update_intocol(&gfiles[i].fltr, gfiles[i].fltr.intocol + ncols);
    ncols = graph_file_usedcols(&gfiles[i]);

    ctx_max_kmers = MAX2(ctx_max_kmers, gfiles[i].num_of_kmers);
  }

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

  // Check graph + paths are compatible
  graphs_paths_compatible(gfiles, num_gfiles, pfiles, num_pfiles);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;
  size_t graph_mem, tmp_path_mem, path_mem;
  char path_mem_str[100];

  // DEV: use threads in memory calculation
  size_t num_of_threads = args->max_work_threads;

  // kmer memory = Edges + paths + 1 bit per colour
  bits_per_kmer = sizeof(Edges)*8 + sizeof(PathIndex)*8 + ncols;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        ctx_max_kmers, false, &graph_mem);

  // Path memory
  tmp_path_mem = path_files_tmp_mem_required(pfiles, num_pfiles);
  path_mem = path_max_mem + tmp_path_mem;

  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  //
  // Open output file
  //
  gzFile gzout;

  if(!args->output_file_set || strcmp("-", args->output_file) == 0)
    gzout = gzdopen(fileno(stdout), "w");
  else
    gzout = gzopen(args->output_file, "w");

  if(gzout == NULL) die("Cannot open output file: %s", args->output_file);

  //
  // Set up memory
  //
  size_t kmer_size = gfiles[0].hdr.kmer_size;

  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, ncols, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(uint8_t));

  // In colour
  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);
  db_graph.node_in_cols = calloc2(bytes_per_col*ncols, sizeof(uint8_t));

  // Paths
  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(PathIndex));
  memset(db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(PathIndex));

  path_store_alloc(&db_graph.pdata, path_max_mem, tmp_path_mem, path_max_usedcols);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++)
    graph_load(&gfiles[i], gprefs, &stats);

  hash_table_print_stats(&db_graph.ht);

  //
  // Load path files
  //
  if(num_pfiles > 0) {
    paths_format_merge(pfiles, num_pfiles, false, &db_graph);
    path_store_reclaim_tmp(&db_graph.pdata);
  }

  //
  // Load reference sequence into a read buffer
  //
  ReadBuffer rbuf;
  readbuf_alloc(&rbuf, 1024);

  read_t r;
  seq_read_alloc(&r);
  for(i = 0; i < num_seq_files; i++) {
    while(seq_read(seq_files[i], &r) > 0) {
      readbuf_add(&rbuf, r);
      seq_read_alloc(&r);
    }
  }
  seq_read_dealloc(&r);

  // Call breakpoints
  breakpoints_call(gzout, &db_graph, rbuf.data, rbuf.len, num_of_threads, args);

  // Finished: do clean up
  gzclose(gzout);

  for(i = 0; i < rbuf.len; i++) seq_read_dealloc(&rbuf.data[i]);

  readbuf_dealloc(&rbuf);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free(db_graph.kmer_paths);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
