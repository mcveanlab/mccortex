#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "graphs_load.h"
#include "graph_writer.h"
#include "clean_graph.h"
#include "db_unitig.h" // for saving length histogram

const char clean_usage[] =
"usage: "CMD" clean [options] <in.ctx> [in2.ctx ...]\n"
"\n"
"  Clean a cortex graph. Joins graphs first, if multiple inputs given.\n"
"  If output graph file is not specified just saves output statistics.\n"
"  If given a multisample graph, cleans each sample against the merged population.\n"
"\n"
"  -h, --help               This help message\n"
"  -q, --quiet              Silence status output normally printed to STDERR\n"
"  -f, --force              Overwrite output files\n"
"  -o, --out <out.ctx>      Save output graph file [required]\n"
"  -m, --memory <mem>       Memory to use\n"
"  -n, --nkmers <kmers>     Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>        Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -N, --ncols <N>          Number of graph colours to use\n"
"  -S, --sort               Output a graph file ordered by kmer\n"
"\n"
"  Cleaning:\n"
"  -T[L], --tips[=L]        Clip tips shorter than <L> kmers [default: auto]\n"
"  -U[X], --unitigs[=X]     Remove low coverage unitigs with median cov < X [default: auto]\n"
"  -B, --fallback <T>       Fall back threshold if we can't pick\n"
"\n"
"  Statistics:\n"
"  -c, --covg-before <out.csv> Save kmer coverage histogram before cleaning\n"
"  -C, --covg-after <out.csv>  Save kmer coverage histogram after cleaning\n"
"  -l, --len-before <out.csv>  Save unitig length histogram before cleaning\n"
"  -L, --len-after <out.csv>   Save unitig length histogram after cleaning\n"
"\n"
"  --unitigs without a threshold, causes a calculated threshold to be used\n"
"  Default: --tips 2*kmer_size --unitigs\n"
"  Set thresholds to zero to turn-off cleaning\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"ncols",        required_argument, NULL, 'N'},
  {"sort",         no_argument,       NULL, 'S'},
// command specific
  {"tips",         optional_argument, NULL, 'T'},
  {"unitigs",      optional_argument, NULL, 'U'},
  {"fallback",     required_argument, NULL, 'B'},
// output
  {"len-before",   required_argument, NULL, 'l'},
  {"len-after",    required_argument, NULL, 'L'},
  {"covg-before",  required_argument, NULL, 'c'},
  {"covg-after",   required_argument, NULL, 'C'},
  {NULL, 0, NULL, 0}
};

// Returns number of kmers in the hash table
static size_t ctx_cleaning_memory(struct MemArgs memargs, bool use_mem_limit,
                                  uint64_t ctx_max_kmers, uint64_t ctx_sum_kmers,
                                  size_t file_ncols, size_t graph_ncols,
                                  bool sort_kmers, size_t *graph_mem_ptr)
{
  bool all_colours_loaded = (file_ncols <= graph_ncols);

  size_t bits_per_kmer, per_col_bits;
  size_t extra_edge_bits, sort_kmers_bits;

  per_col_bits = (sizeof(Covg)+sizeof(Edges)) * 8;
  // We need to store pop edges + sample edges if we haven't loaded all sample
  extra_edge_bits = (all_colours_loaded ? 0 : sizeof(Edges) * 8);
  sort_kmers_bits = (sort_kmers ? sizeof(hkey_t)*8 : 0);

  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  per_col_bits * graph_ncols +
                  extra_edge_bits +
                  sort_kmers_bits;

  return cmd_get_kmers_in_hash(memargs.mem_to_use,
                               memargs.mem_to_use_set,
                               memargs.num_kmers,
                               memargs.num_kmers_set,
                               bits_per_kmer,
                               ctx_max_kmers, ctx_sum_kmers,
                               use_mem_limit, graph_mem_ptr);
}

// Returns number of kmers in the hash table
static size_t ctx_max_cols(struct MemArgs memargs, uint64_t ctx_max_kmers,
                           size_t file_ncols, bool sort_kmers)
{
  // Maximise the number of colours we load to fill the mem
  size_t bits_per_kmer, kmers_in_hash;
  size_t per_col_bits, extra_edge_bits, sort_kmers_bits, ncols;

  kmers_in_hash = ctx_max_kmers / IDEAL_OCCUPANCY;
  per_col_bits = (sizeof(Covg)+sizeof(Edges)) * 8;
  // We need to store pop edges + sample edges if we haven't loaded all sample
  extra_edge_bits = sizeof(Edges) * 8;
  sort_kmers_bits = (sort_kmers ? sizeof(hkey_t)*8 : 0);

  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  per_col_bits * file_ncols +
                  sort_kmers_bits;

  if((bits_per_kmer*kmers_in_hash)/8 <= memargs.mem_to_use) {
    return file_ncols;
  }

  // can't load all colours at once, need extra edge bits
  // remove colour specific bits for now
  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  extra_edge_bits +
                  sort_kmers_bits;

  ncols = (memargs.mem_to_use*8 - bits_per_kmer*kmers_in_hash) /
          (per_col_bits*kmers_in_hash);
  ncols = MAX2(ncols, 1);

  return MIN2(ncols, file_ncols);
}

int ctx_clean(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_ctx_path = NULL;
  bool sort_kmers = false;
  int min_keep_tip = -1, unitig_min = -1; // <0 => default, 0 => noclean
  bool unitig_cleaning = false, tip_cleaning = false;
  uint32_t fallback_thresh = 0;
  const char *len_before_path = NULL, *len_after_path = NULL;
  const char *covg_before_path = NULL, *covg_after_path = NULL;

  // User specified ncols, input colours, how many colours choose to use
  size_t user_ncols = 0, file_ncols = 0, using_ncols = 0;

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o':
        if(out_ctx_path != NULL) cmd_print_usage(NULL);
        out_ctx_path = optarg;
        break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'N': user_ncols = cmd_uint32_nonzero(cmd, optarg); break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'T':
        cmd_check(min_keep_tip<0, cmd);
        min_keep_tip = (optarg != NULL ? (int)cmd_uint32(cmd, optarg) : -1);
        tip_cleaning = true;
        break;
      case 'S': cmd_check(!sort_kmers,cmd); sort_kmers = true; break;
      case 'U':
        cmd_check(unitig_min<0, cmd);
        unitig_min = (optarg != NULL ? cmd_uint32(cmd, optarg) : -1);
        unitig_cleaning = true;
        break;
      case 'B': cmd_check(!fallback_thresh, cmd); fallback_thresh = cmd_uint32_nonzero(cmd, optarg); break;
      case 'l': cmd_check(!len_before_path, cmd); len_before_path = optarg; break;
      case 'L': cmd_check(!len_after_path, cmd); len_after_path = optarg; break;
      case 'c': cmd_check(!covg_before_path, cmd); covg_before_path = optarg; break;
      case 'C': cmd_check(!covg_after_path, cmd); covg_after_path = optarg; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" clean -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;

  if(optind >= argc) cmd_print_usage("Please give input graph files");

  bool doing_cleaning = (unitig_cleaning || tip_cleaning);

  // set default cleaning
  if(!doing_cleaning && out_ctx_path != NULL) {
    unitig_cleaning = tip_cleaning = true;
    doing_cleaning = true;
  }

  // If you ever want to estimate cleaning threshold without outputting
  // a graph, change this to a warning
  if(doing_cleaning && out_ctx_path == NULL) {
    cmd_print_usage("Please specify --out <out.ctx> for cleaned graph");
    // warn("No cleaning being done: you did not specify --out <out.ctx>");
  }

  if(!doing_cleaning && (covg_after_path || len_after_path)) {
    warn("You gave --len-after <out> / --covg-after <out> without "
         "any cleaning (set -U, --unitigs or -t, --tips)");
  }

  if(doing_cleaning && strcmp(out_ctx_path,"-") != 0 &&
     !futil_get_force() && futil_file_exists(out_ctx_path))
  {
    cmd_print_usage("Output file already exists: %s", out_ctx_path);
  }

  if(fallback_thresh && !unitig_cleaning)
    warn("-B, --fallback <T> without --unitigs");

  // Use remaining args as graph files
  char **gfile_paths = argv + optind;
  size_t i, j, num_gfiles = (size_t)(argc - optind);

  // Open graph files
  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t col, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  file_ncols = graph_files_open(gfile_paths, gfiles, num_gfiles,
                                &ctx_max_kmers, &ctx_sum_kmers);

  size_t kmer_size = gfiles[0].hdr.kmer_size;

  // Flatten if we don't have to remember colours / output a graph
  if(out_ctx_path == NULL)
  {
    file_ncols = 1;
    for(i = 0; i < num_gfiles; i++)
      file_filter_flatten(&gfiles[i].fltr, 0);
  }

  if(file_ncols < user_ncols) {
    warn("I only need %zu colour%s ('--ncols %zu' ignored)",
         file_ncols, util_plural_str(file_ncols), user_ncols);
    user_ncols = file_ncols;
  }

  // If no arguments given we default to removing tips < 2*kmer_size
  if(min_keep_tip < 0)
    min_keep_tip = 2 * kmer_size;

  // Warn if any graph files already cleaned
  size_t fromcol;
  ErrorCleaning *cleaning;

  for(i = 0; i < num_gfiles; i++) {
    for(j = 0; j < file_filter_num(&gfiles[i].fltr); j++) {
      fromcol = file_filter_fromcol(&gfiles[i].fltr, j);
      cleaning = &gfiles[i].hdr.ginfo[fromcol].cleaning;
      if(cleaning->cleaned_unitigs && unitig_cleaning) {
        warn("%s:%zu already has unitig cleaning with threshold: <%zu",
             file_filter_path(&gfiles[i].fltr), fromcol,
             (size_t)cleaning->clean_unitigs_thresh);
      }
      if(cleaning->cleaned_tips && tip_cleaning) {
        warn("%s:%zu already has had tip cleaned",
             file_filter_path(&gfiles[i].fltr), fromcol);
      }
    }
  }

  // Print steps
  size_t step = 0;
  status("Actions:\n");
  if(covg_before_path != NULL)
    status("%zu. Saving kmer coverage distribution to: %s", step++, covg_before_path);
  if(len_before_path != NULL)
    status("%zu. Saving unitig length distribution to: %s", step++, len_before_path);
  if(tip_cleaning)
    status("%zu. Cleaning tips shorter than %i nodes", step++, min_keep_tip);
  if(unitig_cleaning) {
    if(unitig_min > 0)
      status("%zu. Cleaning unitigs with coverage < %i", step++, unitig_min);
    if(unitig_min < 0)
      status("%zu. Cleaning unitigs with auto-detected threshold", step++);
  }
  if(covg_after_path != NULL)
    status("%zu. Saving kmer coverage distribution to: %s", step++, covg_after_path);
  if(len_after_path != NULL)
    status("%zu. Saving unitig length distribution to: %s", step++, len_after_path);

  //
  // Decide memory usage
  //
  bool use_mem_limit = (memargs.mem_to_use_set && num_gfiles > 1) || !ctx_max_kmers;
  size_t kmers_in_hash = 0, graph_mem = 0;
  bool all_colours_loaded;

  if(user_ncols)
    using_ncols = user_ncols;
  else
    using_ncols = ctx_max_cols(memargs, ctx_max_kmers, file_ncols, sort_kmers);

  all_colours_loaded = (using_ncols == file_ncols);
  kmers_in_hash = ctx_cleaning_memory(memargs, use_mem_limit,
                                      ctx_max_kmers, ctx_sum_kmers,
                                      file_ncols, using_ncols,
                                      sort_kmers, &graph_mem);

  char max_kmers_str[100];
  ulong_to_str(ctx_max_kmers, max_kmers_str);
  status("[cleaning] %zu input graph%s, max kmers: %s, using %zu colour%s",
         num_gfiles, util_plural_str(num_gfiles), max_kmers_str,
         using_ncols, util_plural_str(using_ncols));

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Check output files are writable
  //
  futil_create_output(out_ctx_path);

  // Does nothing if arg is NULL
  futil_create_output(covg_before_path);
  futil_create_output(covg_after_path);
  futil_create_output(len_before_path);
  futil_create_output(len_after_path);

  // Create db_graph
  // Load as many colours as possible
  // Use an extra set of edge to take intersections
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, using_ncols, using_ncols,
                 kmers_in_hash, DBG_ALLOC_EDGES | DBG_ALLOC_COVGS);

  // Extra edges required to hold union of kept edges
  Edges *edges_union = NULL;
  if(!all_colours_loaded)
    edges_union = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));

  // Load graph into a single colour
  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);

  // Construct cleaned graph header
  GraphFileHeader outhdr;
  memset(&outhdr, 0, sizeof(GraphFileHeader));
  for(i = 0; i < num_gfiles; i++)
    graph_file_merge_header(&outhdr, &gfiles[i]);

  if(!all_colours_loaded)
  {
    db_graph.num_of_cols = db_graph.num_edge_cols = 1;
    SWAP(edges_union, db_graph.col_edges);
    graphs_load_files_flat(gfiles, num_gfiles, gprefs, NULL);
    SWAP(edges_union, db_graph.col_edges);
    db_graph.num_of_cols = db_graph.num_edge_cols = using_ncols;
  }
  else {
    for(i = 0; i < num_gfiles; i++)
      graph_load(&gfiles[i], gprefs, NULL);
  }

  char num_kmers_str[100];
  ulong_to_str(hash_table_nkmers(&db_graph.ht), num_kmers_str);
  status("[cleaning] Total kmers loaded: %s\n", num_kmers_str);

  size_t initial_nkmers = hash_table_nkmers(&db_graph.ht);
  hash_table_print_stats(&db_graph.ht);

  uint8_t *visited = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);
  uint8_t *keep = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

  // Always estimate cleaning threshold
  // if(unitig_min <= 0 || covg_before_path || len_before_path)
  // {
    // Get coverage distribution and estimate cleaning threshold
    int est_min_covg = cleaning_get_threshold(nthreads,
                                              covg_before_path,
                                              len_before_path,
                                              visited, &db_graph);

    if(est_min_covg < 0) status("Cannot find recommended cleaning threshold");
    else status("Recommended cleaning threshold is: %i", est_min_covg);

    // Use estimated threshold if threshold not set
    if(unitig_min < 0) {
      if(fallback_thresh > 0 && est_min_covg < (int)fallback_thresh) {
        status("Using fallback threshold: %i", fallback_thresh);
        unitig_min = fallback_thresh;
      }
      else if(est_min_covg >= 0) unitig_min = est_min_covg;
    }
  // }

  // Die if we failed to find suitable cleaning threshold
  if(unitig_min < 0)
    die("Need cleaning threshold (--unitigs=<D> or --fallback <D>)");

  // Cleaning parameters should now be set (>0) or turned off (==0)
  ctx_assert(unitig_min >= 0);
  ctx_assert(min_keep_tip >= 0);

  if(unitig_cleaning || tip_cleaning)
  {
    // Clean graph of tips (if min_keep_tip > 0) and unitigs (if threshold > 0)
    clean_graph(nthreads, unitig_min, min_keep_tip,
                covg_after_path, len_after_path,
                visited, keep, &db_graph);
  }

  ctx_free(visited);
  ctx_free(keep);

  if(out_ctx_path != NULL)
  {
    // Set output header ginfo cleaned
    for(col = 0; col < using_ncols; col++)
    {
      cleaning = &outhdr.ginfo[col].cleaning;
      cleaning->cleaned_unitigs |= unitig_cleaning;
      cleaning->cleaned_tips |= tip_cleaning;

      // if(tip_cleaning) {
      //   strbuf_append_str(&outhdr.ginfo[col].sample_name, ".tipclean");
      // }

      if(unitig_cleaning) {
        size_t thresh = cleaning->clean_unitigs_thresh;
        thresh = cleaning->cleaned_unitigs ? MAX2(thresh, (uint32_t)unitig_min)
                                          : (uint32_t)unitig_min;
        cleaning->clean_unitigs_thresh = thresh;

        // char name_append[200];
        // sprintf(name_append, ".supclean%zu", thresh);
        // strbuf_append_str(&outhdr.ginfo[col].sample_name, name_append);
      }
    }

    // Print stats on removed kmers
    size_t removed_nkmers = initial_nkmers - hash_table_nkmers(&db_graph.ht);
    double removed_pct = (100.0 * removed_nkmers) / initial_nkmers;
    char removed_str[100], init_str[100];
    ulong_to_str(removed_nkmers, removed_str);
    ulong_to_str(initial_nkmers, init_str);
    status("Removed %s of %s (%.2f%%) kmers", removed_str, init_str, removed_pct);

    // kmers_loaded=true
    graph_writer_merge(out_ctx_path, gfiles, num_gfiles,
                       true, all_colours_loaded,
                       edges_union, &outhdr,
                       sort_kmers, &db_graph);
  }

  ctx_check(hash_table_nkmers(&db_graph.ht) == hash_table_count_kmers(&db_graph.ht));

  // TODO: report kmer coverage for each sample

  graph_header_dealloc(&outhdr);

  for(i = 0; i < num_gfiles; i++) graph_file_close(&gfiles[i]);
  ctx_free(gfiles);

  ctx_free(edges_union);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}