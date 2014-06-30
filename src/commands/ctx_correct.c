#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
#include "correct_reads.h"
#include "read_thread_cmd.h"

const char correct_usage[] =
"usage: "CMD" correct [options] <input.ctx>\n"
"\n"
"  Correct reads against a (population) graph. Uses paths if specified.\n"
"  Bases are printed in lower case if they cannot be corrected.\n"
"\n"
"  -h, --help                 This help message\n"
"  -m, --memory <mem>         Memory to use (e.g. 1M, 20GB)\n"
"  -n, --nkmers <N>           Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>          Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>       Load path file (can specify multiple times)\n"
// Non default:
"  -c, --colour <in:out>      Correct reads from file (supports sam,bam,fq,*.gz\n"
"  -1, --seq <in:out>         Correct reads from file (supports sam,bam,fq,*.gz\n"
"  -2, --seq2 <in1:in2:out>   Correct paired end sequences (output: <out>.{1,2}.fa.gz)\n"
"  -i, --seqi <in.bam:out>    Correct PE reads from a single file\n"
"  -f,--FR -F,--FF            Mate pair orientation [default: FR]\n"
"    -r,--RF -R--RR\n"
"  -w, --oneway               Use one-way gap filling (conservative)\n"
"  -W, --twoway               Use two-way gap filling (liberal)\n"
"  -Q, --fq-threshold <Q>     Filter quality scores [default: 0 (off)]\n"
"  -q, --fq-offset <N>        FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"  -H, --cut-hp <bp>          Breaks reads at homopolymers >= <bp> [default: off]\n"
"  -g, --min-ins <ins>        Minimum insert size for --seq2 [default:"QUOTE_VALUE(DEFAULT_CRTALN_MIN_INS)"]\n"
"  -G, --max-ins <ins>        Maximum insert size for --seq2 [default:"QUOTE_VALUE(DEFAULT_CRTALN_MAX_INS)"]\n"
"\n"
"  -d, --gap-diff-const       -d, -D set parameters for allowable gap lengths\n"
"  -D, --gap-diff-coeff        (gap_exp*D - d) <= gap_actual <= (gap_exp*D + d)\n"
"  -X, --max-context          Number of kmers to use either side of a gap\n"
"  -e, --end-check            Extra check after bridging gap [default: on]\n"
"  -E, --no-end-check         Skip extra check after gap bridging\n"
//
// "  -S, --seq-gaps <out.csv>   Save size distribution of seq gaps bridged\n"
// "  -M, --mp-gaps <out.csv>    Save size distribution of mate pair gaps bridged\n"
//
"\n"
" --seq outputs <out>.fa.gz, --seq2 outputs <out>.1.fa.gz, <out>.2.fa.gz\n"
" --seq must come AFTER two/oneway options. Output may be slightly shuffled.\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",          no_argument,       NULL, 'h'},
  {"out",           required_argument, NULL, 'o'},
  {"memory",        required_argument, NULL, 'm'},
  {"nkmers",        required_argument, NULL, 'n'},
  {"threads",       required_argument, NULL, 't'},
  {"paths",         required_argument, NULL, 'p'},
// command specific
  {"seq",           required_argument, NULL, '1'},
  {"seq2",          required_argument, NULL, '2'},
  {"seqi",          required_argument, NULL, 'i'},
  {"FR",            no_argument,       NULL, 'f'},
  {"FF",            no_argument,       NULL, 'F'},
  {"RF",            no_argument,       NULL, 'r'},
  {"RR",            no_argument,       NULL, 'R'},
  {"oneway",        no_argument,       NULL, 'w'},
  {"twoway",        no_argument,       NULL, 'W'},
  {"fq-cutoff",     required_argument, NULL, 'Q'},
  {"fq-offset",     required_argument, NULL, 'q'},
  {"cut-hp",        required_argument, NULL, 'H'},
  {"min-ins",       no_argument,       NULL, 'g'},
  {"max-ins",       no_argument,       NULL, 'G'},
  {"colour",        required_argument, NULL, 'c'}, // allow --{col,color,colour}
  {"color",         required_argument, NULL, 'c'},
  {"col",           required_argument, NULL, 'c'},
//
  {"gap-diff-const",required_argument, NULL, 'd'},
  {"gap-diff-coeff",required_argument, NULL, 'D'},
  {"max-context",   required_argument, NULL, 'X'},
  {"end-check",     no_argument,       NULL, 'e'},
  {"no-end-check",  no_argument,       NULL, 'E'},
//
  // {"seq-gaps",     required_argument, NULL, 'S'},
  // {"mp-gaps",      required_argument, NULL, 'M'},
  {NULL, 0, NULL, 0}
};


int ctx_correct(int argc, char **argv)
{
  size_t i, j;
  struct ReadThreadCmdArgs args = READ_THREAD_CMD_ARGS_INIT;
  read_thread_args_alloc(&args);
  read_thread_args_parse(&args, argc, argv, longopts, true);

  GraphFileReader *gfile = &args.gfile;
  GPathFileBuffer *gpfiles = &args.gpfiles;
  CorrectAlnInputBuffer *inputs = &args.inputs;
  size_t ctx_total_cols = gfile->hdr.num_of_cols;
  size_t ctx_num_kmers = gfile->num_of_kmers;

  if(args.colour > ctx_total_cols)
    cmd_print_usage("-c %zu is too big [> %zu]", args.colour, ctx_total_cols);

  size_t ctp_usedcols = 0;
  for(i = 0; i < gpfiles->len; i++) {
    if(!file_filter_iscolloaded(&gpfiles->data[i].fltr, args.colour)) {
      cmd_print_usage("Path file doesn't load into colour %zu: %s",
                      args.colour, gpfiles->data[i].fltr.orig_path.buff);
    }
    ctp_usedcols = MAX2(ctp_usedcols, file_filter_usedcols(&gpfiles->data[i].fltr));
  }

  // Check for compatibility between graph files and path files
  graphs_gpaths_compatible(gfile, 1, gpfiles->data, gpfiles->len);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  // 1 bit needed per kmer if we need to keep track of noreseed
  bits_per_kmer = sizeof(GPath)*8 + sizeof(Edges)*8 +
                  (gpfiles->len > 0 ? sizeof(GPath*)*8 : 0) +
                  ctx_num_kmers;

  kmers_in_hash = cmd_get_kmers_in_hash(args.memargs.mem_to_use,
                                        args.memargs.mem_to_use_set,
                                        args.memargs.num_kmers,
                                        args.memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_num_kmers, ctx_num_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = args.memargs.mem_to_use - MIN2(args.memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles->data, gpfiles->len,
                                  ctx_total_cols, rem_mem, false);
  cmd_print_mem(path_mem, "paths");

  // Total memory
  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args.memargs.mem_to_use, total_mem);

  //
  // Check we can read all output files
  //
  // Open output files
  SeqOutput *outputs = ctx_calloc(inputs->len, sizeof(SeqOutput));
  bool output_files_exist = false;

  for(i = 0; i < inputs->len; i++)
  {
    CorrectAlnInput *input = &inputs->data[i];
    input->crt_params.ctxcol = input->crt_params.ctpcol = args.colour;
    SeqOutput *output = &outputs[i];
    seq_output_alloc(output);
    seq_output_set_paths(output, input->out_base,
                         async_task_pe_output(&input->files));
    input->output = output;
    // output check prints warnings and returns true if errors
    output_files_exist |= seq_output_files_exist_check(output);
  }

  // Abandon if some of the output files already exist
  if(output_files_exist) die("Output files already exist");

  // Attempt to open all files
  for(i = 0; i < inputs->len && seq_output_open(&outputs[i]); i++) {}

  // Check if something went wrong - if so remove all output files
  if(i < inputs->len) {
    for(j = 0; j < i; j++) seq_output_delete(&outputs[i]);
    die("Couldn't open output files");
  }

  //
  // Allocate memory
  //

  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile->hdr.kmer_size, ctx_total_cols, 1, kmers_in_hash);

  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);

  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.node_in_cols = ctx_calloc(bytes_per_col * ctx_total_cols, 1);

  // Create a path store that does not tracks path counts
  gpath_reader_alloc_gpstore(gpfiles->data, gpfiles->len, path_mem, false, &db_graph);

  //
  // Load Graph and Path files
  //
  LoadingStats gstats = LOAD_STATS_INIT_MACRO;
  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  // Load graph, print stats, close file
  graph_load(gfile, gprefs, &gstats);
  hash_table_print_stats_brief(&db_graph.ht);
  graph_file_close(gfile);

  // Load path files
  for(i = 0; i < gpfiles->len; i++) {
    gpath_reader_load(&gpfiles->data[i], true, &db_graph);
    gpath_reader_close(&gpfiles->data[i]);
  }

  //
  // Run alignment
  //
  correct_reads(args.num_of_threads, MAX_IO_THREADS,
                inputs->data, inputs->len,
                &db_graph);

  // Close and free output files
  for(i = 0; i < inputs->len; i++) seq_output_dealloc(&outputs[i]);
  ctx_free(outputs);

  read_thread_args_dealloc(&args);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
