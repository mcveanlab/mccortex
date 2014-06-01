#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "loading_stats.h"
#include "graph_info.h"
#include "graph_format.h"
#include "path_format.h"
#include "path_file_filter.h"
#include "graph_file_filter.h"
#include "generate_paths.h"
#include "graph_paths.h"

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
  size_t num_of_threads = DEFAULT_NTHREADS;
  bool mem_to_use_set = false, nkmers_set = false;
  size_t mem_to_use = DEFAULT_MEM, nkmers = DEFAULT_NKMERS;
  char *out_ctp_path = NULL;
  size_t i, max_gap_limit = 0;
  bool use_new_paths = false, clean_paths = false;
  int tmp_thresh, clean_threshold = 0; // 0 => no calling, -1 => auto
  CorrectAlnTask task = {.fq_cutoff = 0, .hp_cutoff = 0,
                         .matedir = READPAIR_FR,
                         .crt_params = CORRECT_PARAMS_DEFAULT,
                         .ptr = NULL};
  uint8_t fq_offset = 0;
  char *dump_seq_sizes = NULL, *dump_mp_sizes = NULL;
  size_t dump_seq_n = 0, dump_mp_n = 0; // how many times are -g -G specified
  PathFileReader tmp_pfile;

  CorrectAlnTaskBuffer tasks;
  correct_aln_task_buf_alloc(&tasks, 16);

  PathFileBuffer pfiles;
  pfile_buf_alloc(&pfiles, 16);

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  long_opts_to_short(longopts, shortopts);
  int used = 1, c;

  printf("shortops: %s\n", shortopts);

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    get_long_opt(longopts, c, cmd);
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': if(out_ctp_path){cmd_print_usage(NULL);} out_ctp_path = optarg; break;
      case 'p':
        tmp_pfile = INIT_PATH_READER;
        path_file_open(&tmp_pfile, optarg, true);
        pfile_buf_add(&pfiles, tmp_pfile);
        break;
      case 't': num_of_threads = cmd_parse_arg_uint32_nonzero(cmd, optarg); break;
      case 'n':
        if(nkmers_set) die("%s specifed more than once", cmd);
        nkmers = cmd_parse_arg_uint32_nonzero(cmd, optarg);
        nkmers_set = true;
        break;
      case 'm':
        if(mem_to_use_set) die("%s specifed more than once", cmd);
        mem_to_use = cmd_parse_arg_mem(cmd, optarg);
        mem_to_use_set = true;
        break;
      case '1':
      case '2':
      case 'i':
        used = 1;
        correct_aln_task_buf_add(&tasks, task);
        asyncio_task_parse(&tasks.data[tasks.len-1].files, c, optarg, fq_offset);
        break;
      case 'f': task.matedir = READPAIR_FR; used = 0; break;
      case 'F': task.matedir = READPAIR_FF; used = 0; break;
      case 'r': task.matedir = READPAIR_RF; used = 0; break;
      case 'R': task.matedir = READPAIR_RR; used = 0; break;
      case 'w': task.crt_params.one_way_gap_traverse = true; used = 0; break;
      case 'W': task.crt_params.one_way_gap_traverse = false; used = 0; break;
      case 'q': fq_offset = cmd_parse_arg_uint8(cmd, optarg); used = 0; break;
      case 'Q': task.fq_cutoff = cmd_parse_arg_uint8(cmd, optarg); used = 0; break;
      case 'H': task.hp_cutoff = cmd_parse_arg_uint8(cmd, optarg); used = 0; break;
      case 'e': task.crt_params.use_end_check = true; used = 0; break;
      case 'E': task.crt_params.use_end_check = false; used = 0; break;
      case 'g': task.crt_params.ins_gap_min = cmd_parse_arg_uint32(cmd, optarg); used = 0; break;
      case 'G': task.crt_params.ins_gap_max = cmd_parse_arg_uint32(cmd, optarg); used = 0; break;
      case 'S': dump_seq_sizes = optarg; dump_seq_n++; break;
      case 'M': dump_mp_sizes = optarg; dump_mp_n++; break;
      case 'u': use_new_paths = true; break;
      case 'c':
        if(optarg == NULL || strcmp(optarg,"auto")) clean_threshold = -1;
        else if(parse_entire_int(optarg,&tmp_thresh) && tmp_thresh >= -1) {
          if(tmp_thresh != -1 && tmp_thresh < 2)
            warn("Ignoring --clean %u (too small < 2)", tmp_thresh);
          else if(tmp_thresh > 255)
            warn("Ignoring --clean %u (too big > 255)", tmp_thresh);
          else
            clean_threshold = tmp_thresh;
        }
        else die("Bad argument for %s <auto|N> where N > 1", cmd);
        clean_paths = (clean_threshold != 0);
        break;
      case 'X': gen_paths_print_contigs = true; break;
      case 'Y': gen_paths_print_paths = true; break;
      case 'Z': gen_paths_print_reads = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" thread -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Check that optind+1 == argc
  if(optind+1 > argc)
    cmd_print_usage("Expected exactly one graph file");
  else if(optind+1 < argc)
    cmd_print_usage("Expected only one graph file. What is this: '%s'", argv[optind]);

  char *graph_path = argv[optind];
  status("Reading graph: %s", graph_path);

  if(!used) die("Ignored arguments after last --seq");

  if(dump_seq_n > 1) die("Cannot specify --seq-gaps <out> more than once");
  if(dump_mp_n > 1) die("Cannot specify --mp-gaps <out> more than once");

  // Check ins_gap_min < ins_gap_max

  for(i = 0; i < tasks.len; i++) {
    CorrectAlnTask *t = &tasks.data[i];
    t->files.ptr = t;
    if(t->crt_params.
      ins_gap_min > t->crt_params.ins_gap_max) {
      die("--min-ins %u is greater than --max-ins %u",
          t->crt_params.ins_gap_min, t->crt_params.ins_gap_max);
    }
    correct_reads_input_print(&tasks.data[i]);
    max_gap_limit = MAX2(max_gap_limit, t->crt_params.ins_gap_max);
  }

  //
  // Open graph graph files
  //
  GraphFileReader gfile = INIT_GRAPH_READER_MACRO;
  graph_file_open(&gfile, graph_path, true);
  file_filter_update_intocol(&gfile.fltr, 0);
  if(graph_file_usedcols(&gfile) > 1)
    die("Please specify a single colour e.g. %s:0", gfile.fltr.file_path.buff);
  size_t ctx_max_kmers = gfile.num_of_kmers;

  //
  // Open path files
  //
  for(i = 0; i < pfiles.len; i++) {
    file_filter_update_intocol(&pfiles.data[i].fltr, 0);
    if(path_file_usedcols(&pfiles.data[i]) > 1) {
      die("Please specify a single colour e.g. %s:0",
          pfiles.data[i].fltr.file_path.buff);
    }
  }

  // Check for compatibility between graph files and path files
  graphs_paths_compatible(&gfile, 1, pfiles.data, pfiles.len);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, total_mem;
  size_t path_mem_req, path_mem;

  bits_per_kmer = sizeof(Edges)*8 +
                  sizeof(PathIndex)*8*(pfiles.len > 0 ? 2 : 1) +
                  2*num_of_threads + // Have traversed
                  1 + // node in colour
                  1; // path store kmer lock

  // false -> don't use mem_to_use to decide how many kmers to store in hash
  // since we need some of that memory for storing paths
  kmers_in_hash = cmd_get_kmers_in_hash2(mem_to_use, mem_to_use_set,
                                         nkmers, nkmers_set,
                                         bits_per_kmer,
                                         ctx_max_kmers, ctx_max_kmers,
                                         false, &graph_mem);

  path_mem_req = path_files_mem_required(pfiles.data, pfiles.len, false, false, 1);
  path_mem = MAX2(mem_to_use - graph_mem, path_mem_req);
  cmd_print_mem(path_mem, "paths");

  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(mem_to_use, total_mem);

  //
  // Open output file
  //
  if(dump_seq_sizes && !futil_is_file_writable(dump_seq_sizes))
    die("Cannot write to file: %s", dump_seq_sizes);
  if(dump_mp_sizes && !futil_is_file_writable(dump_mp_sizes))
    die("Cannot write to file: %s", dump_mp_sizes);

  if(strcmp(out_ctp_path,"-") != 0 && futil_file_exists(out_ctp_path))
    die("Output file already exists: %s", out_ctp_path);

  FILE *fout = !strcmp(out_ctp_path,"-") ? stdout : fopen(out_ctp_path, "w");
  if(fout == NULL) die("Unable to open paths file to write: %s", out_ctp_path);
  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  status("Creating paths file: %s", futil_outpath_str(out_ctp_path));

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, 1, 1, kmers_in_hash);
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
  if(pfiles.len > 0) {
    // Paths loaded into empty colours will update the sample names
    // and add kmers needed
    paths_format_merge(pfiles.data, pfiles.len, true, false,
                       num_of_threads, &db_graph);
  }

  // if no-pickup flags
  if(!use_new_paths) {
    if(pfiles.len > 0) {
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
  pheader.kmer_size = gfile.hdr.kmer_size;
  const char *sample_name = gfile.hdr.ginfo[gfile.fltr.cols[0]].sample_name.buff;
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
                              .empty_colours = true};

  // Load graph
  graph_load(&gfile, gprefs, &gstats);
  hash_table_print_stats_brief(&db_graph.ht);
  graph_file_close(&gfile);

  //
  // Start up the threads, do the work
  //
  GenPathWorker *workers;
  workers = gen_paths_workers_alloc(num_of_threads, &db_graph, fout);

  // Deal with a set of files at once
  size_t start, end;
  for(start = 0; start < tasks.len; start += num_of_threads)
  {
    // Can have different numbers of inputs vs threads
    end = MIN2(tasks.len, start+num_of_threads);
    generate_paths(tasks.data+start, end-start, workers, num_of_threads);
  }

  // Output statistics
  LoadingStats stats = gen_paths_get_stats(workers);
  CorrectAlnStats gapstats =  gen_paths_get_gapstats(workers);

  // Print mp gap size / insert stats to a file
  if(dump_seq_sizes != NULL) {
    correct_aln_stats_dump(dump_seq_sizes,
                             gapstats.gap_err_histgrm, gapstats.histgrm_len,
                             db_graph.kmer_size, false,
                             stats.num_se_reads + stats.num_pe_reads);
  }

  if(stats.num_pe_reads > 0 && dump_mp_sizes != NULL) {
    correct_aln_stats_dump(dump_mp_sizes,
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
  gen_paths_workers_dealloc(workers, num_of_threads);
  path_store_combine_updated_paths(&db_graph.pstore);

  path_store_print_status(&db_graph.pstore);

  // Estimate coverage
  double covg = (double)stats.num_kmers_loaded / db_graph.ht.num_kmers;
  double mean_klen = (double)stats.num_kmers_loaded / stats.contigs_loaded;
  status("Estimate coverage to be %.2f [%zu / %zu], mean len: %.2f",
         covg, stats.num_kmers_loaded, (size_t)db_graph.ht.num_kmers, mean_klen);

  if(clean_paths)
  {
    // Set threshold if not given
    uint8_t threshold;;
    if(clean_threshold == -1) {
      // DEV: calculate threshold
      threshold = 2;
    } else {
      threshold = clean_threshold;
    }

    if(threshold > 1)
      graph_paths_clean(&db_graph, num_of_threads, threshold);
    else
      warn("Path cleaning threshold < 2 has no effect: %i", threshold);
  }

  status("Saving paths to: %s", out_ctp_path);

  // Update header and write
  paths_header_update(&pheader, &db_graph.pstore);
  paths_format_write_header(&pheader, fout);
  paths_format_write_optimised_paths(&db_graph, fout);
  fclose(fout);

  paths_header_dealloc(&pheader);

  for(i = 0; i < tasks.len; i++) asyncio_task_close(&tasks.data[i].files);
  correct_aln_task_buf_dealloc(&tasks);

  ctx_free(db_graph.node_in_cols);
  ctx_free(db_graph.col_edges);
  path_store_dealloc(&db_graph.pstore);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < pfiles.len; i++) path_file_close(&pfiles.data[i]);
  pfile_buf_dealloc(&pfiles);

  return EXIT_SUCCESS;
}
