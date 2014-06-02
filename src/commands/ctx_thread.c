#include "global.h"
#include "commands.h"
#include "file_util.h"
#include "db_graph.h"
#include "loading_stats.h"
#include "graph_format.h"
#include "path_format.h"
#include "path_file_filter.h"
#include "graph_file_filter.h"
#include "generate_paths.h"
#include "graph_paths.h"
#include "read_thread_cmd.h"

const char thread_usage[] =
"usage: "CMD" thread [options] <in.ctx>\n"
"\n"
"  Thread reads through the graph.  Save to file <out.ctp>.  <pop.ctx> can should\n"
"  have everyone in it and can be a pooled graph (with only 1 colour).  Samples\n"
"  are loaded from <in.ctx> files one at a time.\n"
"\n"
"  -h, --help               This help message\n"
"  -o, --out <out.ctp>      Save output file [required]\n"
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
"  -Q, --fq_threshold <Q>   Filter quality scores [default: 0 (off)]\n"
"  -q, --fq_offset <N>      FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"  -H, --cut_hp <bp>        Breaks reads at homopolymers >= <bp> [default: off]\n"
"  -e, --end-check          Extra check after bridging gap [default: on]\n"
"  -E, --no-end-check       Skip extra check after gap bridging\n"
"  -g, --min-ins <ins>      Minimum insert size for --seq2 [default:0]\n"
"  -G, --max-ins <ins>      Maximum insert size for --seq2 [default:"QUOTE_MACRO(DEFAULT_MAX_INS)"]\n"
"  -S, --seq-gaps <out.csv> Save size distribution of seq gaps bridged\n"
"  -M, --mp-gaps <out.csv>  Save size distribution of mate pair gaps bridged\n"
"  -u, --use-new-paths      Use paths as they are being added (higher err rate) [default: no]\n"
"  -c, --clean [auto|N]     Threshold at auto or <N> and remove redundant paths\n"
"\n"
"  Debugging Options: Probably best not to touch these\n"
"    -X,--print-contigs -Y,--print-paths -Z,--print-reads\n"
"\n"
"  When loading existing paths with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into. See `"CMD" pjoin` to combine .ctp files\n"
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
  {"seq-gaps",     required_argument, NULL, 'S'},
  {"mp-gaps",      required_argument, NULL, 'M'},
  {"use-new-paths",optional_argument, NULL, 'u'},
  {"clean",        optional_argument, NULL, 'c'},
// Debug options
  {"print-contigs",no_argument,       NULL, 'X'},
  {"print-paths",  no_argument,       NULL, 'Y'},
  {"print-reads",  no_argument,       NULL, 'Z'},
  {NULL, 0, NULL, 0}
};


int ctx_thread(int argc, char **argv)
{
  struct ReadThreadCmdArgs args = READ_THREAD_CMD_ARGS_INIT;
  read_thread_args_alloc(&args);
  read_thread_args_parse(&args, argc, argv, longopts, false);

  GraphFileReader *gfile = &args.gfile;
  PathFileBuffer *pfiles = &args.pfiles;
  CorrectAlnInputBuffer *inputs = &args.inputs;

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  bits_per_kmer = sizeof(Edges)*8 +
                  sizeof(PathIndex)*8*(pfiles->len > 0 ? 2 : 1) +
                  2*args.num_of_threads + // Have traversed
                  1 + // node in colour
                  1; // path store kmer lock

  // false -> don't use mem_to_use to decide how many kmers to store in hash
  // since we need some of that memory for storing paths
  kmers_in_hash = cmd_get_kmers_in_hash2(args.memargs.mem_to_use,
                                         args.memargs.mem_to_use_set,
                                         args.memargs.num_kmers,
                                         args.memargs.num_kmers_set,
                                         bits_per_kmer,
                                         gfile->num_of_kmers,
                                         gfile->num_of_kmers,
                                         false, &graph_mem);

  path_mem = path_files_mem_required(pfiles->data, pfiles->len, false, false, 0);
  path_mem = MAX2(args.memargs.mem_to_use - graph_mem, path_mem);
  cmd_print_mem(path_mem, "paths");

  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args.memargs.mem_to_use, total_mem);

  //
  // Open output file
  //
  if(args.dump_seq_sizes && !futil_is_file_writable(args.dump_seq_sizes))
    die("Cannot write to file: %s", args.dump_seq_sizes);
  if(args.dump_mp_sizes && !futil_is_file_writable(args.dump_mp_sizes))
    die("Cannot write to file: %s", args.dump_mp_sizes);

  if(strcmp(args.out_ctp_path,"-") != 0 && futil_file_exists(args.out_ctp_path))
    die("Output file already exists: %s", args.out_ctp_path);

  FILE *fout = !strcmp(args.out_ctp_path,"-") ? stdout : fopen(args.out_ctp_path, "w");
  if(fout == NULL) die("Unable to open paths file to write: %s", args.out_ctp_path);
  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  status("Creating paths file: %s", futil_outpath_str(args.out_ctp_path));

  //
  // Allocate memory
  //
  dBGraph db_graph;
  size_t kmer_size = gfile->hdr.kmer_size;
  db_graph_alloc(&db_graph, kmer_size, 1, 1, kmers_in_hash);
  kmers_in_hash = db_graph.ht.capacity;

  // Edges
  db_graph.col_edges = ctx_calloc(kmers_in_hash, sizeof(Edges));

  // Path store
  path_store_alloc(&db_graph.pstore, path_mem, true, kmers_in_hash, 1);

  // Keep counts in an extra byte per path
  db_graph.pstore.extra_bytes = 1;

  // path kmer locks for multithreaded access
  db_graph.pstore.kmer_locks = ctx_calloc(roundup_bits2bytes(kmers_in_hash), 1);

  // Load existing paths
  if(pfiles->len > 0) {
    // Paths loaded into empty colours will update the sample names
    // and add kmers needed
    paths_format_merge(pfiles->data, pfiles->len, true, false,
                       args.num_of_threads, &db_graph);
  }

  // if no-pickup flags
  if(!args.use_new_paths) {
    if(pfiles->len > 0) {
      // Copy current paths over to path set to be updated
      size_t mem = kmers_in_hash * sizeof(PathIndex);
      db_graph.pstore.kmer_paths_read = ctx_malloc(mem);
      memcpy(db_graph.pstore.kmer_paths_read, db_graph.pstore.kmer_paths_write, mem);
    }
    else
      db_graph.pstore.kmer_paths_read = NULL;
  }

  // Set up paths header. This is for the output file we are creating
  PathFileHeader pheader = INIT_PATH_FILE_HDR_MACRO;
  paths_header_alloc(&pheader, 1);

  pheader.num_of_cols = 1;
  pheader.kmer_size = kmer_size;
  const char *sample_name = gfile->hdr.ginfo[gfile->fltr.cols[0]].sample_name.buff;
  strbuf_set(&pheader.sample_names[0], sample_name);

  // 2. reduce number of graph colours
  db_graph_realloc(&db_graph, 1, 1);

  db_graph.node_in_cols = ctx_calloc(roundup_bits2bytes(kmers_in_hash), 1);

  // Setup for loading graphs graph
  LoadingStats gstats;
  loading_stats_init(&gstats);

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  // Load graph, print stats, close file
  graph_load(gfile, gprefs, &gstats);
  hash_table_print_stats_brief(&db_graph.ht);
  graph_file_close(gfile);

  //
  // Start up the threads, do the work
  //
  GenPathWorker *workers;
  workers = gen_paths_workers_alloc(args.num_of_threads, &db_graph, fout);

  // Deal with a set of files at once
  size_t start, end;
  for(start = 0; start < inputs->len; start += args.num_of_threads)
  {
    // Can have different numbers of inputs vs threads
    end = MIN2(inputs->len, start+args.num_of_threads);
    generate_paths(inputs->data+start, end-start, workers, args.num_of_threads);
  }

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
  path_store_combine_updated_paths(&db_graph.pstore);

  path_store_print_status(&db_graph.pstore);

  // Estimate coverage
  double covg = (double)stats.num_kmers_loaded / db_graph.ht.num_kmers;
  double mean_klen = (double)stats.num_kmers_loaded / stats.contigs_loaded;
  status("Estimate coverage to be %.2f [%zu / %zu], mean len: %.2f",
         covg, stats.num_kmers_loaded, (size_t)db_graph.ht.num_kmers, mean_klen);

  if(args.clean_paths)
  {
    // Set threshold if not given
    uint8_t threshold;;
    if(args.clean_threshold == -1) {
      // DEV: calculate threshold
      threshold = 2;
    } else {
      threshold = args.clean_threshold;
    }

    if(threshold > 1)
      graph_paths_clean(&db_graph, args.num_of_threads, threshold);
    else
      warn("Path cleaning threshold < 2 has no effect: %i", threshold);
  }

  status("Saving paths to: %s", args.out_ctp_path);

  // Update header and write
  paths_header_update(&pheader, &db_graph.pstore);
  paths_format_write_header(&pheader, fout);
  paths_format_write_optimised_paths(&db_graph, fout);
  fclose(fout);

  paths_header_dealloc(&pheader);

  read_thread_args_dealloc(&args);

  ctx_free(db_graph.node_in_cols);
  ctx_free(db_graph.col_edges);
  path_store_dealloc(&db_graph.pstore);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
