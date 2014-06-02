#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_paths.h"
#include "correct_reads.h"

// DEV: print read sequence in lower case instead of N

const char correct_usage[] =
"usage: "CMD" correct [options] <input.ctx> [...]\n"
"\n"
"  Correct reads against a (population) graph.\n"
"\n"
"  -h, --help                 This help message\n"
"  -m, --memory <mem>         Memory to use (e.g. 1M, 20GB)\n"
"  -n, --nkmers <N>           Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>          Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>       Load path file (can specify multiple times)\n"
// Non default:
"  -1, --seq <in:out>         Correct reads from file (supports sam,bam,fq,*.gz\n"
"  -2, --seq2 <in1:in2:out>   Correct paired end sequences (output: <out>.{1,2}.fa.gz)\n"
"  -i, --seqi <in.bam:out>    Correct PE reads from a single file\n"
"  -f,--FR -F,--FF            Mate pair orientation [default: FR]\n"
"    -r,--RF -R--RR\n"
"  -w, --oneway               Use one-way gap filling (conservative)\n"
"  -W, --twoway               Use two-way gap filling (liberal)\n"
"  -Q, --fq_threshold <Q>     Filter quality scores [default: 0 (off)]\n"
"  -q, --fq_offset <N>        FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"  -H, --cut_hp <bp>          Breaks reads at homopolymers >= <bp> [default: off]\n"
"  -e, --end-check            Extra check after bridging gap [default: on]\n"
"  -E, --no-end-check         Skip extra check after gap bridging\n"
"  -g, --min-ins <ins>        Minimum insert size for --seq2 [default:0]\n"
"  -G, --max-ins <ins>        Maximum insert size for --seq2 [default:"QUOTE_MACRO(DEFAULT_MAX_INS)"]\n"
// "  -S, --seq-gaps <out.csv>   Save size distribution of seq gaps bridged\n"
// "  -M, --mp-gaps <out.csv>    Save size distribution of mate pair gaps bridged\n"
"\n"
" --seq outputs <out>.fa.gz, --seq2 outputs <out>.1.fa.gz, <out>.2.fa.gz\n"
" --seq must come AFTER two/oneway options. Output may be slightly shuffled.\n"
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
// command specific
  {"seq",          required_argument, NULL, '1'},
  {"seq2",         required_argument, NULL, '2'},
  {"seqi",         required_argument, NULL, 'i'},
  {"FR",           no_argument,       NULL, 'f'},
  {"FF",           no_argument,       NULL, 'F'},
  {"RF",           no_argument,       NULL, 'r'},
  {"RR",           no_argument,       NULL, 'R'},
  {"oneway",       no_argument,       NULL, 'w'},
  {"twoway",       no_argument,       NULL, 'W'},
  {"fq_cutoff",    required_argument, NULL, 'Q'},
  {"fq_offset",    required_argument, NULL, 'q'},
  {"cut_hp",       required_argument, NULL, 'H'},
  {"end-check",    no_argument,       NULL, 'e'},
  {"no-end-check", no_argument,       NULL, 'E'},
  {"min-ins",      no_argument,       NULL, 'g'},
  {"max-ins",      no_argument,       NULL, 'G'},
  // {"seq-gaps",     required_argument, NULL, 'S'},
  // {"mp-gaps",      required_argument, NULL, 'M'},
  {NULL, 0, NULL, 0}
};

int ctx_correct(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that we have at least 2 arguments

  size_t max_ninputs = argc / 2, num_inputs = 0;
  CorrectAlnTask *inputs = ctx_malloc(max_ninputs * sizeof(CorrectAlnTask));
  size_t i, j;
  size_t num_threads = args->max_work_threads;
  size_t max_io_threads = args->max_io_threads;
  int argi; // arg index to continue from

  // Load args
  // argi = load_args(argc, argv, inputs, &num_inputs);
  argi = correct_reads_parse(argc, argv, false, true,
                             inputs, &num_inputs,
                             NULL, NULL);

  if(argi == argc) cmd_print_usage("Expected at least one graph file");
  size_t num_gfiles = argc - argi;
  char **graph_paths = &argv[argi];

  //
  // Open graph gfiles
  //
  GraphFileReader gfiles[num_gfiles];
  size_t ctx_total_cols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ctx_total_cols = graph_files_open(graph_paths, gfiles, num_gfiles,
                                    &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Open path files
  //
  size_t num_pfiles = args->num_ctp_files;
  PathFileReader pfiles[num_pfiles];
  size_t ctp_max_usedcols = 0;

  for(i = 0; i < num_pfiles; i++) {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], args->ctp_files[i], true);
    ctp_max_usedcols = MAX2(ctp_max_usedcols, path_file_usedcols(&pfiles[i]));
  }

  // Check for compatibility between graph files and path files
  graphs_paths_compatible(gfiles, num_gfiles, pfiles, num_pfiles);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  // 1 bit needed per kmer if we need to keep track of noreseed
  bits_per_kmer = sizeof(Edges)*8 + ctx_max_kmers + sizeof(uint64_t)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, &graph_mem);

  // Paths memory
  path_mem = path_files_mem_required(pfiles, num_pfiles, false, false, 0);
  cmd_print_mem(path_mem, "paths");

  // Total memory
  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args->mem_to_use, total_mem);

  //
  // Check we can read all output files
  //
  // Loop over inputs, change input->ptr from char* to CorrectedOutput
  CorrectedOutput outputs[num_inputs];

  for(i = 0; i < num_inputs && corrected_output_open(&outputs[i], &inputs[i]); i++)
    inputs[i].ptr = &outputs[i];

  // Check if something went wrong
  if(i < num_inputs) {
    for(j = 0; j < i; j++)
      corrected_output_delete(&outputs[i]);
    die("Couldn't open output files");
  }

  // futil_is_file_writable() creates the output file - we have to delete these
  // if we then fail
  for(i = 0; i < num_inputs; i++) {
    if(futil_file_exists((char*)inputs[i].ptr)) {
      for(j = 0; j < i; j++) unlink((char*)inputs[j].ptr);
      die("File already exists: %s", (char*)inputs[i].ptr);
    }
    if(!futil_is_file_writable((char*)inputs[i].ptr)) {
      for(j = 0; j < i; j++) unlink((char*)inputs[j].ptr);
      die("Cannot write to output file: %s", (char*)inputs[i].ptr);
    }
  }

  // Allocate
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ctx_total_cols, 1, kmers_in_hash);

  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);

  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.node_in_cols = ctx_calloc(bytes_per_col * ctx_total_cols, 1);

  // Paths
  path_store_alloc(&db_graph.pstore, path_mem, false,
                   db_graph.ht.capacity, ctp_max_usedcols);

  //
  // Load Graph and Path files
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

  // Load path files (does nothing if num_fpiles == 0)
  paths_format_merge(pfiles, num_pfiles, false, false,
                     args->max_work_threads, &db_graph);

  for(i = 0; i < num_pfiles; i++) path_file_close(&pfiles[i]);

  //
  // Run alignment
  //
  correct_reads(num_threads, max_io_threads, inputs, outputs, num_inputs,
                &db_graph);

  ctx_free(inputs);
  ctx_free(db_graph.col_edges);
  ctx_free(db_graph.node_in_cols);

  path_store_dealloc(&db_graph.pstore);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
