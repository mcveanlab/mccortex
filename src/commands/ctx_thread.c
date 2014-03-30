#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "loading_stats.h"
#include "graph_info.h"
#include "graph_format.h"
#include "path_format.h"
#include "generate_paths.h"
#include "graph_paths.h"

const char thread_usage[] =
"usage: "CMD" thread [options] <in.ctx> [in2.ctx[:cols] ...]\n"
"  Thread reads through the graph.  Save to file <out.ctp>.  <pop.ctx> can should\n"
"  have everyone in it and can be a pooled graph (with only 1 colour).  Samples\n"
"  are loaded from <in.ctx> files one at a time.\n"
"\n"
"  Options:\n"
"    -m <mem>                   How much memory to use\n"
"    -n <kmers>                 How many entries in the hash table\n"
"    -p <in.ctp>                Load existing path files first (multiple allowed)\n"
"    --out <out.ctp>            Save output file [required]\n"
"    --col <colour>             Colour to thread through\n"
"    --seq <in.fa>              Thread reads from file (supports sam,bam,fq,*.gz)\n"
"    --seq2 <in.1.fq> <in.2.fq> Thread paired end reads\n"
"    --minIns <ins>             Minimum insert size for --seq2 [default:0]\n"
"    --maxIns <ins>             Maximum insert size for --seq2 [default:"QUOTE_MACRO(DEFAULT_MAX_INS)"]\n"
"    --oneway                   Use one-way gap filling [default]\n"
"    --twoway                   Use two-way gap filling\n"
"    --fq_threshold <fq>        FASTQ quality threshold\n"
"    --fq_offset <qual>         FASTQ quality score offset\n"
"    --cut_hp <N>               Cut reads afer <N> consecutive bases\n"
"    --FR --FF --RF --RR        Mate pair orientation [default: FR] (with --keep_pcr)\n"
"    --seqgaps <out.csv>        Save size distribution of seq gaps bridged\n"
"    --mpgaps <out.csv>         Save size distribution of mate pair gaps bridged\n"
"\n"
"  Debugging Options:\n"
"    --printcontigs --printpaths --printreads   Probably best not to touch these\n"
"\n"
"  When loading existing paths with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into. See `"CMD" pjoin` to combine .ctp files\n";



static void get_binary_and_colour(const GraphFileReader *files, size_t num_gfiles,
                                  size_t col, size_t *file_idx, size_t *col_idx)
{
  size_t i, n = 0;
  for(i = 0; i < num_gfiles; i++) {
    if(n + graph_file_outncols(&files[i]) > col) {
      *col_idx = col - n; *file_idx = i; return;
    }
    n += graph_file_outncols(&files[i]);
  }
  die("Colour is greater than sum of graph colours [%zu > %zu]", col, n);
}

int ctx_thread(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;

  size_t max_tasks = (size_t)argc/2;
  CorrectAlnReadsTask *tasks = malloc2(max_tasks * sizeof(CorrectAlnReadsTask));
  size_t i, j, num_tasks, num_work_threads = args->max_work_threads;
  int argi; // arg index to continue from
  char *dump_seq_sizes = NULL, *dump_mp_sizes = NULL;

  argi = correct_reads_parse(argc, argv, tasks, &num_tasks,
                             true, false, &dump_seq_sizes, &dump_mp_sizes);

  for(i = 0; i < num_tasks; i++) {
    tasks[i].crt_params.ctxcol = 0;
    correct_reads_input_print(&tasks[i]);
  }

  if(argi == argc) cmd_print_usage("Not enough arguments");

  const char *out_ctp_path = args->output_file;
  size_t num_gfiles = (size_t)(argc - argi);
  char **graph_paths = argv + argi;

  //
  // Open graph graph files
  //
  GraphFileReader gfiles[num_gfiles];
  size_t ctx_max_kmers = 0, ctx_total_cols = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    gfiles[i] = INIT_GRAPH_READER;
    graph_file_open(&gfiles[i], graph_paths[i], true);

    file_filter_update_intocol(&gfiles[i].fltr,
                               gfiles[i].fltr.intocol + ctx_total_cols);
    ctx_total_cols = graph_file_usedcols(&gfiles[i]);
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

  // Check for compatibility between graph files and path files
  graphs_paths_compatible(gfiles, num_gfiles, pfiles, num_pfiles);

  // Check for path colours >= number of output path colours
  size_t max_gap_limit = 0;

  for(i = 0; i < num_tasks; i++) {
    max_gap_limit = MAX2(max_gap_limit, tasks[i].crt_params.ins_gap_max);

    if(tasks[i].crt_params.ctpcol >= ctx_total_cols) {
      cmd_print_usage("--col <C> >= number of path colours [%zu >= %zu]",
                      tasks[i].crt_params.ctpcol, ctx_total_cols);
    }
  }

  size_t total_cols = MAX2(ctx_total_cols, path_max_usedcols);
  status("Creating paths file with %zu colour%s: %s",
         total_cols, total_cols != 1 ? "s" : "", out_ctp_path);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, total_mem;
  size_t tmp_path_mem, path_mem_req, path_mem_used, main_path_mem;
  char path_mem_str[100];

  bits_per_kmer = sizeof(Edges)*8 +
                  sizeof(PathIndex)*8*(num_pfiles > 0 ? 2 : 1) +
                  2*num_work_threads + // Have traversed
                  1 + // node in colour
                  1; // path store kmer lock

  // false -> don't use mem_to_use to decide how many kmers to store in hash
  // since we need some of that memory for storing paths
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, ctx_max_kmers,
                                        false, &graph_mem);

  // Path Memory
  tmp_path_mem = path_files_tmp_mem_required(pfiles, num_pfiles);
  path_mem_req = path_max_mem + tmp_path_mem;
  path_mem_used = MAX2(args->mem_to_use - graph_mem, path_mem_req);
  main_path_mem = path_mem_used - tmp_path_mem;

  bytes_to_str(path_mem_used, 1, path_mem_str);
  status("[memory] paths: %s", path_mem_str);

  total_mem = graph_mem + path_mem_used;
  cmd_check_mem_limit(args, total_mem);

  //
  // Open output file
  //
  if(dump_seq_sizes && !futil_is_file_writable(dump_seq_sizes))
    die("Cannot write to file: %s", dump_seq_sizes);
  if(dump_mp_sizes && !futil_is_file_writable(dump_mp_sizes))
    die("Cannot write to file: %s", dump_mp_sizes);

  if(futil_file_exists(out_ctp_path))
    die("Output file already exists: %s", out_ctp_path);

  FILE *fout = fopen(out_ctp_path, "w");
  if(fout == NULL) die("Unable to open paths file to write: %s", out_ctp_path);
  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, total_cols, 1, kmers_in_hash);
  kmers_in_hash = db_graph.ht.capacity;

  // Edges
  db_graph.col_edges = calloc2(kmers_in_hash, sizeof(Edges));

  // Path store
  // use total_cols instead of path_max_usedcols since we are
  // loading then ADDING more paths (which may need new colours)
  path_store_alloc(&db_graph.pstore, main_path_mem, tmp_path_mem,
                   kmers_in_hash, total_cols);

  // path kmer locks for multithreaded access
  db_graph.pstore.kmer_locks = calloc2(roundup_bits2bytes(kmers_in_hash), 1);

  // if no-pickup flags
  if(num_pfiles > 0)
    db_graph.pstore.kmer_paths = malloc2(kmers_in_hash * sizeof(PathIndex));
  else
    db_graph.pstore.kmer_paths = NULL;
  //

  // 1. Merge graph file headers into the graph
  size_t intocol, fromcol;
  for(i = 0; i < num_gfiles; i++) {
    for(j = 0; j < gfiles[i].fltr.ncols; j++) {
      intocol = graph_file_intocol(&gfiles[i], j);
      fromcol = graph_file_fromcol(&gfiles[i], j);
      graph_info_merge(&db_graph.ginfo[intocol],
                       &gfiles[i].hdr.ginfo[fromcol]);
    }
  }

  // Load existing paths
  if(num_pfiles > 0) {
    // Paths loaded into empty colours will update the sample names
    // and add kmers needed
    paths_format_merge(pfiles, num_pfiles, true, &db_graph);
    path_store_reclaim_tmp(&db_graph.pstore);

    // Copy current paths over to path set to be updated
    memcpy(db_graph.pstore.kmer_paths_update, db_graph.pstore.kmer_paths,
           kmers_in_hash * sizeof(PathIndex));
  }
  else {
    // Don't pick up any paths
    db_graph.pstore.kmer_paths = NULL;
  }

  // Set up paths header. This is for the output file we are creating
  PathFileHeader pheader = INIT_PATH_FILE_HDR_MACRO;
  paths_header_alloc(&pheader, total_cols);

  pheader.num_of_cols = (uint32_t)total_cols;
  pheader.kmer_size = gfiles[0].hdr.kmer_size;

  // Set new path header sample names from graph header
  for(i = 0; i < total_cols; i++)
    strbuf_set(&pheader.sample_names[i], db_graph.ginfo[i].sample_name.buff);

  // 2. reduce number of graph colours
  db_graph_realloc(&db_graph, 1, 1);

  db_graph.node_in_cols = calloc2(roundup_bits2bytes(kmers_in_hash), 1);

  // Setup for loading graphs graph
  LoadingStats gstats;
  loading_stats_init(&gstats);

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = false};

  //
  // Start up the threads, do the work
  //
  GenPathWorker *workers;
  workers = gen_paths_workers_alloc(num_work_threads, &db_graph);

  // ... Send jobs ...
  size_t start, end, ctpcol, prev_ctpcol = SIZE_MAX;
  size_t fileidx = 0, colidx = 0;
  for(start = 0; start < num_tasks; start = end, prev_ctpcol = ctpcol)
  {
    // Load graph
    ctpcol = tasks[start].crt_params.ctpcol;

    if(ctpcol != prev_ctpcol)
    {
      if(start > 0) {
        // wipe colour 0: just reset edges
        memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
        // DO NOT EMPTY THE HASH TABLE stupid
        // uncommenting this messes up keeping track of kmer->path matching
        // hash_table_empty(&db_graph.ht);
      }

      get_binary_and_colour(gfiles, num_gfiles, ctpcol, &fileidx, &colidx);
      graph_load_colour(&gfiles[fileidx], gprefs, &gstats, colidx, 0);

      hash_table_print_stats_brief(&db_graph.ht);
    }

    // Get list of input files to read
    end = start+1;
    while(end < num_tasks && end-start < args->max_io_threads &&
          tasks[end].crt_params.ctpcol == ctpcol) end++;

    generate_paths(tasks+start, end-start, workers, num_work_threads);
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
  size_t num_gaps_paths_disagreed = gapstats.num_gaps_disagreed;
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

  PathStore *pstore = &db_graph.pstore;
  size_t num_path_bytes = (size_t)(pstore->next - pstore->store);
  char kmers_str[100], paths_str[100], mem_str[100], col_paths_str[100];
  ulong_to_str(pstore->num_kmers_with_paths, kmers_str);
  ulong_to_str(pstore->num_of_paths, paths_str);
  bytes_to_str(num_path_bytes, 1, mem_str);
  ulong_to_str(pstore->num_col_paths, col_paths_str);

  // ins_gap, err_gap no longer allocated after this line
  gen_paths_workers_dealloc(workers, num_work_threads);
  path_store_combine_updated_paths(&db_graph.pstore);

  status("Saving paths: %s paths, %s path-bytes, %s kmers, coloured paths: %s",
         paths_str, mem_str, kmers_str, col_paths_str);

  // Update header and write
  paths_header_update(&pheader, pstore);
  paths_format_write_header(&pheader, fout);
  paths_format_write_optimised_paths(&db_graph, fout);
  fclose(fout);

  for(i = 0; i < num_tasks; i++) {
    if(tasks[i].file1 != NULL) seq_close(tasks[i].file1);
    if(tasks[i].file2 != NULL) seq_close(tasks[i].file2);
  }

  free(tasks);

  free(db_graph.node_in_cols);
  free(db_graph.col_edges);

  paths_header_dealloc(&pheader);

  path_store_dealloc(&db_graph.pstore);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_gfiles; i++) graph_file_dealloc(&gfiles[i]);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  status("Paths written to: %s", out_ctp_path);

  return EXIT_SUCCESS;
}
