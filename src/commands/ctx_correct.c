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
"  -h, --help               This help message\n"
"  -q, --quiet              Silence status output normally printed to STDERR\n"
"  -f, --force              Overwrite output files\n"
"  -m, --memory <mem>       Memory to use (e.g. 1M, 20GB)\n"
"  -n, --nkmers <N>         Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>        Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>     Load path file (can specify multiple times)\n"
"\n"
"  Input:\n"
"  -1, --seq <in:out>       Correct reads (output: <out>.fa.gz)\n"
"  -2, --seq2 <in1:in2:out> Correct reads (output: <out>{,.1,.2}.fa.gz)\n"
"  -i, --seqi <in.bam:out>  Correct reads (output: <out>{,.1,.2}.fa.gz)\n"
"  -M, --matepair <orient>  Mate pair orientation: FF,FR,RF,RR [default: FR]\n"
"  -Q, --fq-cutoff <Q>      Filter quality scores [default: 0 (off)]\n"
"  -O, --fq-offset <N>      FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"  -H, --cut-hp <bp>        Breaks reads at homopolymers >= <bp> [default: off]\n"
"  -l, --frag-len-min <bp>  Min fragment size for --seq2 [default:"QUOTE_VALUE(DEFAULT_CRTALN_FRAGLEN_MIN)"]\n"
"  -L, --frag-len-max <bp>  Max fragment size for --seq2 [default:"QUOTE_VALUE(DEFAULT_CRTALN_FRAGLEN_MAX)"]\n"
"\n"
"  Path Params:\n"
"  -w, --one-way            Use one-way gap filling (conservative)\n"
"  -W, --two-way            Use two-way gap filling (liberal)\n"
"  -d, --gap-diff-const     -d, -D set parameters for allowable gap lengths\n"
"  -D, --gap-diff-coeff       (gap_exp*D - d) <= gap_actual <= (gap_exp*D + d)\n"
"  -X, --max-context        Number of kmers to use either side of a gap\n"
"  -e, --end-check          Extra check after bridging gap [default: on]\n"
"  -E, --no-end-check       Skip extra check after gap bridging\n"
"  -G, --gap-hist <o.csv>   Save size distribution of sequence gaps bridged\n"
"  -F, --frag-hist <o.csv>  Save size distribution of PE fragments recovered\n"
"\n"
"  -c, --colour <in:out>    Correct reads from file (supports sam,bam,fq,*.gz\n"
"\n"
" --seq outputs <out>.fa.gz, --seq2 outputs <out>.1.fa.gz, <out>.2.fa.gz\n"
" --seq must come AFTER two/oneway options. Output may be slightly shuffled.\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",          no_argument,       NULL, 'h'},
  {"memory",        required_argument, NULL, 'm'},
  {"nkmers",        required_argument, NULL, 'n'},
  {"threads",       required_argument, NULL, 't'},
  {"paths",         required_argument, NULL, 'p'},
  {"force",         no_argument,       NULL, 'f'},
// command specific
  {"seq",           required_argument, NULL, '1'},
  {"seq2",          required_argument, NULL, '2'},
  {"seqi",          required_argument, NULL, 'i'},
  {"matepair",      required_argument, NULL, 'M'},
  {"fq-cutoff",     required_argument, NULL, 'Q'},
  {"fq-offset",     required_argument, NULL, 'O'},
  {"cut-hp",        required_argument, NULL, 'H'},
  {"min-frag-len",  required_argument, NULL, 'l'},
  {"max-frag-len",  required_argument, NULL, 'L'},
//
  {"one-way",       no_argument,       NULL, 'w'},
  {"two-way",       no_argument,       NULL, 'W'},
  {"gap-diff-const",required_argument, NULL, 'd'},
  {"gap-diff-coeff",required_argument, NULL, 'D'},
  {"max-context",   required_argument, NULL, 'X'},
  {"end-check",     no_argument,       NULL, 'e'},
  {"no-end-check",  no_argument,       NULL, 'E'},
  {"gap-hist",      required_argument, NULL, 'G'},
  {"frag-hist",     required_argument, NULL, 'F'},
//
  {"colour",        required_argument, NULL, 'c'}, // allow --{col,color,colour}
  {"color",         required_argument, NULL, 'c'},
  {"col",           required_argument, NULL, 'c'},
//
  // {"print-contigs", no_argument,       NULL, 'x'},
  // {"print-paths",   no_argument,       NULL, 'y'},
  // {"print-reads",   no_argument,       NULL, 'z'},
  {NULL, 0, NULL, 0}
};


int ctx_correct(int argc, char **argv)
{
  size_t i;
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
                  ctx_total_cols; // in colour

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

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  // Total memory
  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args.memargs.mem_to_use, total_mem);

  //
  // Check we can write all output files
  //
  // Open output files
  SeqOutput *outputs = ctx_calloc(inputs->len, sizeof(SeqOutput));
  bool err_occurred = false;

  // Always print with FASTA format
  for(i = 0; i < inputs->len && !err_occurred; i++)
  {
    CorrectAlnInput *input = &inputs->data[i];
    input->crt_params.ctxcol = input->crt_params.ctpcol = args.colour;
    bool is_pe = asyncio_task_is_pe(&input->files);
    err_occurred = !seqout_open(&outputs[i], input->out_base, SEQ_FMT_FASTQ, is_pe);
    input->output = &outputs[i];
  }

  // Abandon if some of the output files already exist
  if(err_occurred) {
    for(i = 0; i < inputs->len; i++)
      seqout_close(&outputs[i], true);
    die("Error creating output files");
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
  correct_reads(inputs->data, inputs->len,
                args.dump_seq_sizes, args.dump_frag_sizes,
                args.nthreads, &db_graph);

  // Close and free output files
  for(i = 0; i < inputs->len; i++)
    seqout_close(&outputs[i], false);
  ctx_free(outputs);

  // Closes input files
  read_thread_args_dealloc(&args);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
