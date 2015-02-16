#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "graph_format.h"
#include "clean_graph.h"
#include "supernode.h" // for saving length histogram

const bool use_supernode_covg = false;

const char clean_usage[] =
"usage: "CMD" clean [options] <in.ctx> [in2.ctx ...]\n"
"\n"
"  Clean a cortex graph. Joins graphs first, if multiple inputs given.\n"
"  If neither -t or -s specified, just saves output statistics.\n"
"\n"
"  -h, --help                  This help message\n"
"  -q, --quiet                 Silence status output normally printed to STDERR\n"
"  -f, --force                 Overwrite output files\n"
"  -o, --out <out.ctx>         Save output graph file [required]\n"
"  -m, --memory <mem>          Memory to use\n"
"  -n, --nkmers <kmers>        Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>           Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -N, --ncols <N>             Number of graph colours to use\n"
"\n"
"  Cleaning:\n"
"  -T, --tips <L>              Clip tips shorter than <L> kmers\n"
"  -S[T], --supernodes[=T]     Remove low coverage supernode with coverage < T [default: auto]\n"
"  -d, --kdepth <C>            kmer depth: (depth*(R-Kmersize+1)/R); R = read length\n"
"\n"
"  Statistics:\n"
"  -c, --covg-before <out.csv> Save kmer coverage histogram before cleaning\n"
"  -C, --covg-after <out.csv>  Save kmer coverage histogram after cleaning\n"
"  -l, --len-before <out.csv>  Save supernode length histogram before cleaning\n"
"  -L, --len-after <out.csv>   Save supernode length histogram after cleaning\n"
"\n"
"  --supernodes without a threshold, causes a calculated threshold to be used\n"
"  Default: --tips 2*kmer_size --supernodes\n"
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
// command specific
  {"tips",         required_argument, NULL, 'T'},
  {"supernodes",   optional_argument, NULL, 'S'},
  {"kdepth",       required_argument, NULL, 'd'},
// output
  {"len-before",   required_argument, NULL, 'l'},
  {"len-after",    required_argument, NULL, 'L'},
  {"covg-before",  required_argument, NULL, 'c'},
  {"covg-after",   required_argument, NULL, 'C'},
  {NULL, 0, NULL, 0}
};

int ctx_clean(int argc, char **argv)
{
  size_t nthreads = 0, use_ncols = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_ctx_path = NULL;
  bool tip_cleaning = false, supernode_cleaning = false;
  size_t min_keep_tip = 0;
  Covg threshold = 0;
  double seq_depth = 0;
  const char *len_before_path = NULL, *len_after_path = NULL;
  const char *covg_before_path = NULL, *covg_after_path = NULL;

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
      case 'N': use_ncols = cmd_uint32_nonzero(cmd, optarg); break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'T':
        cmd_check(!tip_cleaning, cmd);
        min_keep_tip = cmd_uint32_nonzero(cmd, optarg);
        tip_cleaning = true;
        break;
      case 'S':
        cmd_check(!supernode_cleaning, cmd);
        if(optarg != NULL) threshold = cmd_uint32_nonzero(cmd, optarg);
        supernode_cleaning = true;
        break;
      case 'd': cmd_check(seq_depth <= 0, cmd); seq_depth = cmd_udouble_nonzero(cmd, optarg); break;
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

  // Default behaviour
  if(!tip_cleaning && !supernode_cleaning) {
    if(out_ctx_path != NULL)
      supernode_cleaning = tip_cleaning = true; // do both
    else
      warn("No cleaning being done: you did not specify --out <out.ctx>");
  }

  bool doing_cleaning = (supernode_cleaning || tip_cleaning);

  if(doing_cleaning && out_ctx_path == NULL) {
    cmd_print_usage("Please specify --out <out.ctx> for cleaned graph");
  }

  if(!supernode_cleaning && seq_depth > 0)
    cmd_print_usage("--kdepth <C> not needed if not cleaning with --supernodes");

  if(supernode_cleaning && threshold != 0 && seq_depth > 0)
    cmd_print_usage("--supernodes <T> arg makes --kdepth <D> redundant");

  if(!doing_cleaning && (covg_after_path || len_after_path)) {
    cmd_print_usage("You gave --len-after <out> / --covg-after <out> without "
                    "any cleaning (set -s, --supernodes or -t, --tips)");
  }

  if(doing_cleaning && strcmp(out_ctx_path,"-") != 0 &&
     !futil_get_force() && futil_file_exists(out_ctx_path))
  {
    cmd_print_usage("Output file already exists: %s", out_ctx_path);
  }

  // Use remaining args as graph files
  char **gfile_paths = argv + optind;
  size_t i, j, num_gfiles = (size_t)(argc - optind);

  // Open graph files
  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(gfile_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  size_t kmer_size = gfiles[0].hdr.kmer_size;

  // default to one colour for now
  if(use_ncols == 0) use_ncols = 1;

  // Flatten if we don't have to remember colours / output a graph
  if(!doing_cleaning)
  {
    ncols = use_ncols = 1;
    for(i = 0; i < num_gfiles; i++)
      file_filter_flatten(&gfiles[i].fltr, 0);
  }

  if(ncols < use_ncols) {
    warn("I only need %zu colour%s ('--ncols %zu' ignored)",
         ncols, util_plural_str(ncols), use_ncols);
    use_ncols = ncols;
  }

  char max_kmers_str[100];
  ulong_to_str(ctx_max_kmers, max_kmers_str);
  status("%zu input graph%s, max kmers: %s, using %zu colours",
         num_gfiles, util_plural_str(num_gfiles), max_kmers_str, use_ncols);

  // If no arguments given we default to removing tips < 2*kmer_size
  if(tip_cleaning && min_keep_tip == 0)
    min_keep_tip = 2 * kmer_size;

  // Warn if any graph files already cleaned
  size_t fromcol, intocol;
  ErrorCleaning *cleaning;

  for(i = 0; i < num_gfiles; i++) {
    for(j = 0; j < file_filter_num(&gfiles[i].fltr); j++) {
      fromcol = file_filter_fromcol(&gfiles[i].fltr, j);
      cleaning = &gfiles[i].hdr.ginfo[fromcol].cleaning;
      if(cleaning->cleaned_snodes && supernode_cleaning) {
        warn("%s:%zu already has supernode cleaning with threshold: <%zu",
             file_filter_path(&gfiles[i].fltr), fromcol,
             (size_t)cleaning->clean_snodes_thresh);
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
    status("%zu. Saving supernode length distribution to: %s", step++, len_before_path);
  if(tip_cleaning)
    status("%zu. Cleaning tips shorter than %zu nodes", step++, min_keep_tip);
  if(supernode_cleaning && threshold > 0)
    status("%zu. Cleaning supernodes with coverage < %u", step++, threshold);
  if(supernode_cleaning && threshold == 0)
    status("%zu. Cleaning supernodes with auto-detected threshold", step++);
  if(covg_after_path != NULL)
    status("%zu. Saving kmer coverage distribution to: %s", step++, covg_after_path);
  if(len_after_path != NULL)
    status("%zu. Saving supernode length distribution to: %s", step++, len_after_path);

  //
  // Decide memory usage
  //
  bool all_colours_loaded = (ncols <= use_ncols);
  bool use_mem_limit = (memargs.mem_to_use_set && num_gfiles > 1) || !ctx_max_kmers;

  size_t kmers_in_hash, bits_per_kmer, graph_mem;
  size_t per_kmer_per_col_bits = (sizeof(BinaryKmer)+sizeof(Covg)+sizeof(Edges)) * 8;
  size_t pop_edges_per_kmer_bits = (all_colours_loaded ? 0 : sizeof(Edges) * 8);

  bits_per_kmer = per_kmer_per_col_bits * use_ncols + pop_edges_per_kmer_bits;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        use_mem_limit, &graph_mem);

  // Maximise the number of colours we load to fill the mem
  size_t max_usencols = (memargs.mem_to_use*8 - pop_edges_per_kmer_bits * kmers_in_hash) /
                        (per_kmer_per_col_bits * kmers_in_hash);
  use_ncols = MIN2(max_usencols, ncols);

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
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, use_ncols, use_ncols,
                 kmers_in_hash, DBG_ALLOC_COVGS);

  // Edges is a special case
  size_t num_edges = db_graph.ht.capacity * (use_ncols + !all_colours_loaded);
  db_graph.col_edges = ctx_calloc(num_edges, sizeof(Edges));

  // Load graph into a single colour
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = false};

  // Construct cleaned graph header
  GraphFileHeader outhdr;
  memset(&outhdr, 0, sizeof(GraphFileHeader));
  outhdr.version = CTX_GRAPH_FILEFORMAT;
  outhdr.kmer_size = db_graph.kmer_size;
  outhdr.num_of_cols = ncols;
  outhdr.num_of_bitfields = (db_graph.kmer_size*2+63)/64;
  graph_header_alloc(&outhdr, ncols);

  // Merge info into header
  size_t gcol = 0;
  for(i = 0; i < num_gfiles; i++) {
    for(j = 0; j < file_filter_num(&gfiles[i].fltr); j++, gcol++) {
      fromcol = file_filter_fromcol(&gfiles[i].fltr, j);
      intocol = file_filter_intocol(&gfiles[i].fltr, j);
      graph_info_merge(&outhdr.ginfo[intocol], &gfiles[i].hdr.ginfo[fromcol]);
    }
  }

  if(ncols > use_ncols) {
    graph_files_load_flat(gfiles, num_gfiles, gprefs, &stats);
  } else {
    for(i = 0; i < num_gfiles; i++)
      graph_load(&gfiles[i], gprefs, &stats);
  }

  char num_kmers_str[100];
  ulong_to_str(db_graph.ht.num_kmers, num_kmers_str);
  status("Total kmers loaded: %s\n", num_kmers_str);

  size_t initial_nkmers = db_graph.ht.num_kmers;
  hash_table_print_stats(&db_graph.ht);

  uint8_t *visited = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);
  uint8_t *keep = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

  if(threshold == 0 || covg_before_path || len_before_path) {
    // Get coverage distribution and estimate cleaning threshold
    int est_threshold = cleaning_get_threshold(nthreads, use_supernode_covg,
                                               seq_depth,
                                               covg_before_path, len_before_path,
                                               visited, &db_graph);

    // Use estimated threshold if threshold not set
    if(threshold == 0) {
      // Die if we failed to find suitable cleaning threshold
      if(est_threshold < 0) die("Cannot find suitable cleaning threshold");
      else threshold = est_threshold;
    }
  }

  if(doing_cleaning) {
    // Clean graph of tips (if min_keep_tip > 0) and supernodes (if threshold > 0)
    clean_graph(nthreads, use_supernode_covg, threshold, min_keep_tip,
                covg_after_path, len_after_path,
                visited, keep, &db_graph);
  }

  ctx_free(visited);
  ctx_free(keep);

  if(doing_cleaning)
  {
    // Output graph file
    Edges *intersect_edges = NULL;
    bool kmers_loaded = true;
    size_t col, thresh;

    // Set output header ginfo cleaned
    for(col = 0; col < ncols; col++)
    {
      cleaning = &outhdr.ginfo[col].cleaning;
      cleaning->cleaned_snodes |= supernode_cleaning;
      cleaning->cleaned_tips |= tip_cleaning;

      if(tip_cleaning) {
        strbuf_append_str(&outhdr.ginfo[col].sample_name, ".tipclean");
      }

      if(supernode_cleaning) {
        thresh = cleaning->clean_snodes_thresh;
        thresh = cleaning->cleaned_snodes ? MAX2(thresh, threshold) : threshold;
        cleaning->clean_snodes_thresh = thresh;

        char name_append[200];
        sprintf(name_append, ".supclean%zu", thresh);
        strbuf_append_str(&outhdr.ginfo[col].sample_name, name_append);
      }
    }

    if(!all_colours_loaded)
    {
      // We haven't loaded all the colours
      // intersect_edges are edges to mask with
      // resets graph edges
      intersect_edges = db_graph.col_edges;
      db_graph.col_edges += db_graph.ht.capacity;
    }

    // Print stats on removed kmers
    size_t removed_nkmers = initial_nkmers - db_graph.ht.num_kmers;
    double removed_pct = (100.0 * removed_nkmers) / initial_nkmers;
    char removed_str[100], init_str[100];
    ulong_to_str(removed_nkmers, removed_str);
    ulong_to_str(initial_nkmers, init_str);
    status("Removed %s of %s (%.2f%%) kmers", removed_str, init_str, removed_pct);

    graph_files_merge(out_ctx_path, gfiles, num_gfiles,
                      kmers_loaded, all_colours_loaded,
                      intersect_edges, &outhdr, &db_graph);

    // Swap back
    if(!all_colours_loaded)
      db_graph.col_edges = intersect_edges;
  }

  ctx_check(db_graph.ht.num_kmers == hash_table_count_kmers(&db_graph.ht));

  graph_header_dealloc(&outhdr);

  for(i = 0; i < num_gfiles; i++) graph_file_close(&gfiles[i]);
  ctx_free(gfiles);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
