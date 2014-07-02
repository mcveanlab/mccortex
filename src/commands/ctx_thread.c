#include "global.h"
#include "commands.h"
#include "file_util.h"
#include "db_graph.h"
#include "loading_stats.h"
#include "graph_format.h"
#include "graph_file_reader.h"
#include "generate_paths.h"
#include "read_thread_cmd.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
#include "gpath_save.h"

const char thread_usage[] =
"usage: "CMD" thread [options] <in.ctx>\n"
"\n"
"  Thread reads through the graph.  Save to file <out.ctp>.  <pop.ctx> can should\n"
"  have everyone in it and can be a pooled graph (with only 1 colour).  Samples\n"
"  are loaded from <in.ctx> files one at a time.\n"
"\n"
"  -h, --help               This help message\n"
"  -o, --out <out.ctp.gz>   Save output file [required]\n"
"  -m, --memory <mem>       Memory to use (e.g. 1M, 20GB)\n"
"  -n, --nkmers <N>         Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>        Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>     Load path file (can specify multiple times)\n"
// Non default:
"  -1, --seq <in.fa>        Thread reads from file (supports sam,bam,fq,*.gz\n"
"  -2, --seq2 <in1:in2>     Thread paired end sequences\n"
"  -i, --seqi <in.bam>      Thread PE reads from a single file\n"
"  -f,--FR -F,--FF          Mate pair orientation [default: FR]\n"
"    -r,--RF -R--RR\n"
"  -w, --oneway             Use one-way gap filling (conservative)\n"
"  -W, --twoway             Use two-way gap filling (liberal)\n"
"  -Q, --fq-threshold <Q>   Filter quality scores [default: 0 (off)]\n"
"  -q, --fq-offset <N>      FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"  -H, --cut-hp <bp>        Breaks reads at homopolymers >= <bp> [default: off]\n"
"  -g, --min-ins <ins>      Minimum insert size for --seq2 [default:"QUOTE_VALUE(DEFAULT_CRTALN_MIN_INS)"]\n"
"  -G, --max-ins <ins>      Maximum insert size for --seq2 [default:"QUOTE_VALUE(DEFAULT_CRTALN_MAX_INS)"]\n"
"  -S, --seq-gaps <out.csv> Save size distribution of seq gaps bridged\n"
"  -M, --mp-gaps <out.csv>  Save size distribution of mate pair gaps bridged\n"
"  -u, --use-new-paths      Use paths as they are being added (higher err rate) [default: no]\n"
"\n"
"  -a, --aln-const            -a, -A set parameters for allowable gap lengths\n"
"  -A, --aln-coeff              (gap_exp*A - a) <= gap_actual <= (gap_exp*A + a)\n"
"  -X, --max-context          Number of kmers to use either side of a gap\n"
"  -e, --end-check            Extra check after bridging gap [default: on]\n"
"  -E, --no-end-check         Skip extra check after gap bridging\n"
"\n"
"  Debugging Options: Probably best not to touch these\n"
"    -x,--print-contigs -y,--print-paths -z,--print-reads\n"
"\n"
"  When loading existing paths with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into. See `"CMD" pjoin` to combine .ctp files\n"
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
  {"min-ins",       required_argument, NULL, 'g'},
  {"max-ins",       required_argument, NULL, 'G'},
  {"seq-gaps",      required_argument, NULL, 'S'},
  {"mp-gaps",       required_argument, NULL, 'M'},
  {"use-new-paths", required_argument, NULL, 'u'},
//
  {"gap-diff-const",required_argument, NULL, 'd'},
  {"gap-diff-coeff",required_argument, NULL, 'D'},
  {"max-context",   required_argument, NULL, 'X'},
  {"end-check",     no_argument,       NULL, 'e'},
  {"no-end-check",  no_argument,       NULL, 'E'},
// Debug options
  {"print-contigs", no_argument,       NULL, 'x'},
  {"print-paths",   no_argument,       NULL, 'y'},
  {"print-reads",   no_argument,       NULL, 'z'},
  {NULL, 0, NULL, 0}
};


int ctx_thread(int argc, char **argv)
{
  struct ReadThreadCmdArgs args = READ_THREAD_CMD_ARGS_INIT;
  read_thread_args_alloc(&args);
  read_thread_args_parse(&args, argc, argv, longopts, false);

  GraphFileReader *gfile = &args.gfile;
  GPathFileBuffer *gpfiles = &args.gpfiles;
  CorrectAlnInputBuffer *inputs = &args.inputs;
  size_t i;

  //
  // Decide on memory
  //
  size_t ncols = 1;
  size_t bits_per_kmer, kmers_in_hash, graph_mem, total_mem;
  size_t path_hash_mem, path_store_mem, path_mem;
  bool sep_path_list = (!args.use_new_paths && gpfiles->len > 0);

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 + sizeof(GPath*)*8 +
                  2 * args.num_of_threads; // Have traversed

  // false -> don't use mem_to_use to decide how many kmers to store in hash
  // since we need some of that memory for storing paths
  kmers_in_hash = cmd_get_kmers_in_hash(args.memargs.mem_to_use,
                                        args.memargs.mem_to_use_set,
                                        args.memargs.num_kmers,
                                        args.memargs.num_kmers_set,
                                        bits_per_kmer,
                                        gfile->num_of_kmers,
                                        gfile->num_of_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t min_path_mem = 0, max_path_mem = 0;
  gpath_reader_max_mem_req(gpfiles->data, gpfiles->len,
                           ncols, kmers_in_hash,
                           true, sep_path_list, true,
                           &min_path_mem, &max_path_mem);

  // Maximise path memory
  path_mem = min_path_mem;
  if(graph_mem + path_mem < args.memargs.mem_to_use)
    path_mem = args.memargs.mem_to_use - graph_mem;

  path_hash_mem = path_mem / 3;
  cmd_print_mem(path_hash_mem, "paths hash");
  path_store_mem = (2*path_mem) / 3;
  cmd_print_mem(path_store_mem, "paths store");

  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args.memargs.mem_to_use, total_mem);

  //
  // Open output file
  //
  if(args.dump_seq_sizes && !futil_is_file_writable(args.dump_seq_sizes))
    die("Cannot write to file: %s", args.dump_seq_sizes);
  if(args.dump_mp_sizes && !futil_is_file_writable(args.dump_mp_sizes))
    die("Cannot write to file: %s", args.dump_mp_sizes);

  gzFile gzout = futil_gzopen_output(args.out_ctp_path);

  status("Creating paths file: %s", futil_outpath_str(args.out_ctp_path));

  //
  // Allocate memory
  //
  dBGraph db_graph;
  size_t kmer_size = gfile->hdr.kmer_size;
  db_graph_alloc(&db_graph, kmer_size, ncols, ncols, kmers_in_hash);
  kmers_in_hash = db_graph.ht.capacity;

  // Edges
  db_graph.col_edges = ctx_calloc(kmers_in_hash, sizeof(Edges));

  // Split path memory 2:1 between store and hash
  // Create a path store that tracks path counts
  gpath_store_alloc(&db_graph.gpstore,
                    db_graph.num_of_cols, db_graph.ht.capacity,
                    0, path_store_mem, true, sep_path_list);

  // Create path hash table for fast lookup
  gpath_hash_alloc(&db_graph.gphash, &db_graph.gpstore, path_hash_mem);

  // Load existing paths
  // Paths loaded into empty colours will update the sample names
  // and add kmers needed
  for(i = 0; i < gpfiles->len; i++)
    gpath_reader_load(&gpfiles->data[i], false, &db_graph);

  if(args.use_new_paths) {
    status("Using paths as they are added (risky)");
  } else {
    gpath_store_split_read_write(&db_graph.gpstore);
    status("Not using new paths as they are added (safe)");
  }

  db_graph.node_in_cols = ctx_calloc(roundup_bits2bytes(kmers_in_hash), 1);

  // Setup for loading graphs graph
  LoadingStats gstats;
  loading_stats_init(&gstats);

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = false}; // already loaded paths

  // Load graph, print stats, close file
  graph_load(gfile, gprefs, &gstats);
  hash_table_print_stats_brief(&db_graph.ht);
  graph_file_close(gfile);

  //
  // Start up the threads, do the work
  //
  GenPathWorker *workers;
  workers = gen_paths_workers_alloc(args.num_of_threads, &db_graph);

  // Deal with a set of files at once
  size_t start, end;
  for(start = 0; start < inputs->len; start += MAX_IO_THREADS)
  {
    // Can have different numbers of inputs vs threads
    end = MIN2(inputs->len, start+MAX_IO_THREADS);
    generate_paths(inputs->data+start, end-start, workers, args.num_of_threads);
  }

  // Print memory statistics
  gpath_hash_print_stats(&db_graph.gphash);
  gpath_store_print_stats(&db_graph.gpstore);

  // Output statistics
  LoadingStats stats = gen_paths_get_stats(workers);
  CorrectAlnStats gapstats =  gen_paths_get_gapstats(workers);

  // Print mp gap size / insert stats to a file
  if(args.dump_seq_sizes != NULL) {
    correct_aln_stats_dump(args.dump_seq_sizes,
                             gapstats.gap_err_histgrm, gapstats.histgrm_len,
                             db_graph.kmer_size, false,
                             stats.num_se_reads + stats.num_pe_reads);
  }

  if(stats.num_pe_reads > 0 && args.dump_mp_sizes != NULL) {
    correct_aln_stats_dump(args.dump_mp_sizes,
                           gapstats.gap_ins_histgrm, gapstats.histgrm_len,
                           db_graph.kmer_size, true, stats.num_pe_reads);
  }

  // Path Stats
  size_t num_gap_attempts = gapstats.num_gap_attempts;
  size_t num_gap_successes = gapstats.num_gap_successes;
  size_t num_gaps_paths_disagreed = gapstats.num_paths_disagreed;
  size_t num_gaps_too_short = gapstats.num_gaps_too_short;
  char num_gap_attempts_str[100], num_gap_successes_str[100];
  char num_gaps_paths_disagree_str[100], num_gaps_too_short_str[100];
  ulong_to_str(num_gap_attempts, num_gap_attempts_str);
  ulong_to_str(num_gap_successes, num_gap_successes_str);
  ulong_to_str(num_gaps_paths_disagreed, num_gaps_paths_disagree_str);
  ulong_to_str(num_gaps_too_short, num_gaps_too_short_str);

  status("[gaps] traversals succeeded: %s / %s (%.2f%%)",
         num_gap_successes_str, num_gap_attempts_str,
         (100.0 * num_gap_successes) / num_gap_attempts);
  status("[gaps] failed path check: %s / %s (%.2f%%)",
         num_gaps_paths_disagree_str, num_gap_attempts_str,
         (100.0 * num_gaps_paths_disagreed) / num_gap_attempts);
  status("[gaps] too short: %s / %s (%.2f%%)",
         num_gaps_too_short_str, num_gap_attempts_str,
         (100.0 * num_gaps_too_short) / num_gap_attempts);

  char se_num_str[100], pe_num_str[100], sepe_num_str[100];
  ulong_to_str(stats.num_se_reads, se_num_str);
  ulong_to_str(stats.num_pe_reads / 2, pe_num_str);
  ulong_to_str(stats.num_se_reads + stats.num_pe_reads, sepe_num_str);
  status("[stats] single reads: %s; read pairs: %s; total: %s",
         se_num_str, pe_num_str, sepe_num_str);

  // ins_gap, err_gap no longer allocated after this line
  gen_paths_workers_dealloc(workers, args.num_of_threads);

  // Estimate coverage
  double covg = (double)stats.num_kmers_loaded / db_graph.ht.num_kmers;
  double mean_klen = (double)stats.num_kmers_loaded / stats.contigs_loaded;
  status("Estimate coverage to be %.2f [%zu / %zu], mean len: %.2f",
         covg, stats.num_kmers_loaded, (size_t)db_graph.ht.num_kmers, mean_klen);

  status("Saving paths to: %s", args.out_ctp_path);

  cJSON *hdrs[gpfiles->len];
  for(i = 0; i < gpfiles->len; i++) hdrs[i] = gpfiles->data[i].json;

  // Write output file
  gpath_save(gzout, args.out_ctp_path, hdrs, gpfiles->len, &db_graph);
  gzclose(gzout);

  gpath_checks_all_paths(&db_graph);
  read_thread_args_dealloc(&args);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
