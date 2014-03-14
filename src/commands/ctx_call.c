#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_paths.h"
#include "graph_walker.h"
#include "bubble_caller.h"
#include "db_node.h"

// Long flanks help us map calls
// increasing allele length can be costly
#define DEFAULT_MAX_FLANK 1000
#define DEFAULT_MAX_ALLELE 300

const char call_usage[] =
"usage: "CMD" call [options] <in.ctx> [in2.ctx ...]\n"
"  Find bubbles in the graph, which are potential variants.\n"
"\n"
"  Options:\n"
"    -m <memory> | -t <threads> | -p <paths.ctp>\n"
"    -p <in.ctp>          Load path file (can specify multiple times)\n"
"    -o <out.bubbles.gz>  Output file [required]\n"
"    --haploid <col>      Colour is haploid, can use repeatedly [e.g. ref colour]\n"
"    --maxallele <len>    Max bubble branch length in kmers [default: "QUOTE_MACRO(DEFAULT_MAX_ALLELE)"]\n"
"    --maxflank <len>     Max flank length in kmers [default: "QUOTE_MACRO(DEFAULT_MAX_FLANK)"]\n"
"\n"
"  When loading path files with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into.\n";

int ctx_call(CmdArgs *args)
{
  int argi, argc = args->argc;
  char **argv = args->argv;

  const char *out_path = args->output_file;
  size_t num_of_threads = args->max_work_threads;
  size_t i, haploid_cols[argc], num_haploid = 0;
  size_t max_allele_len = DEFAULT_MAX_ALLELE, max_flank_len = DEFAULT_MAX_FLANK;

  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++)
  {
    if(strcmp(argv[argi],"--haploid") == 0) {
      if(argi + 1 == argc ||
         !parse_entire_size(argv[argi+1], &haploid_cols[num_haploid])) {
        cmd_print_usage("--haploid <col> requires an int arg");
      }
      num_haploid++; argi++;
    }
    else if(strcmp(argv[argi],"--maxallele") == 0) {
      if(argi+1 == argc || !parse_entire_size(argv[argi+1], &max_allele_len) ||
         max_allele_len == 0)  {
        cmd_print_usage("--maxallele <col> requires an +ve integer argument");
      }
      argi++;
    }
    else if(strcmp(argv[argi],"--maxflank") == 0) {
      if(argi+1 == argc || !parse_entire_size(argv[argi+1], &max_flank_len) ||
         max_flank_len == 0)  {
        cmd_print_usage("--maxflank <col> requires an +ve integer argument");
      }
      argi++;
    }
    else {
      cmd_print_usage("Unknown arg: %s", argv[argi]);
    }
  }

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
  // Check haploid colours are valid
  //
  for(i = 0; i < num_haploid; i++) {
    if(haploid_cols[i] >= ncols) {
      cmd_print_usage("--haploid <col> is greater than max colour [%zu > %zu]",
                      haploid_cols[i], ncols-1);
    }
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t tmp_path_mem, path_mem, thread_mem;
  char path_mem_str[100], thread_mem_str[100];

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +
  // visitedfw/rv(2bits/thread)

  bits_per_kmer = sizeof(Edges)*8 + sizeof(PathIndex)*8 +
                  ncols + 2*num_of_threads;

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        ctx_max_kmers, false, &graph_mem);

  // Thread memory
  thread_mem = roundup_bits2bytes(kmers_in_hash) * 2;
  bytes_to_str(thread_mem * num_of_threads, 1, thread_mem_str);
  status("[memory] (of which threads: %zu x %zu = %s)\n",
          num_of_threads, thread_mem, thread_mem_str);

  // Path Memory
  tmp_path_mem = path_files_tmp_mem_required(pfiles, num_pfiles);
  path_mem = path_max_mem + tmp_path_mem;

  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  size_t total_mem = graph_mem + thread_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  //
  // Open output file
  //
  gzFile gzout;

  if(strcmp("-", out_path) == 0)
    gzout = gzdopen(fileno(stdout), "w");
  else
    gzout = gzopen(out_path, "w");

  if(gzout == NULL) die("Cannot open output file: %s", out_path);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash);

  // Edges merged into one colour
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(uint8_t));

  // In colour
  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);
  db_graph.node_in_cols = calloc2(bytes_per_col*ncols, sizeof(uint8_t));

  // Paths
  path_store_alloc(&db_graph.pstore, path_max_mem, tmp_path_mem,
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
    gprefs.empty_colours = false;
  }

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  if(num_pfiles > 0) {
    paths_format_merge(pfiles, num_pfiles, false, &db_graph);
    path_store_reclaim_tmp(&db_graph.pstore);
  }

  // Now call variants
  bubble_caller_print_header(&db_graph, gzout, out_path, args);
  invoke_bubble_caller(&db_graph, gzout, num_of_threads,
                       max_allele_len, max_flank_len, haploid_cols, num_haploid);

  status("  saved to: %s\n", out_path);
  gzclose(gzout);

  free(db_graph.col_edges);
  free(db_graph.node_in_cols);

  path_store_dealloc(&db_graph.pstore);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_gfiles; i++) graph_file_dealloc(&gfiles[i]);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  return EXIT_SUCCESS;
}
