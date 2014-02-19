#include "global.h"

#include "tools.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "graph_format.h"
#include "clean_graph.h"

const char clean_usage[] =
"usage: "CMD" clean [options] <in.ctx> [in2.ctx ...]\n"
" Clean a cortex binary. Joins binaries first, if multiple given\n"
" Clips tips before doing supernode thresholding (when doing both [default]).\n"
"\n"
" Options:\n"
"  --memory <mem>        Memory to use\n"
"  --nkmers <hash-size>  Kmers in the hash table (e.g. 1G ~ 1 billion)\n"
"  --ncols <colour>      Number of samples in memory at once (speedup)\n"
"  --tips <L>            Clip tips shorter than <L> kmers\n"
"  --supernodes          Remove low coverage supernode. Additional options:\n"
"    --kdepth <C>        kmer depth: (depth*(R-Kmersize+1)/R); R = read length\n"
"    --threshold <T>     Cleaning threshold, remove supnodes where [coverage < T]\n"
"  --covgs <out.csv>     Dump covg distribution before cleaning to a CSV file\n"
"  --out <out.ctx>       Save output file\n"
"\n"
" Default: --tips 2*k-1 --supernodes\n";

// size of coverage histogram 2^11 = 2048
#define HISTSIZE 2048

#ifdef CTXVERBOSE
  #define DEBUG_SUPERNODE 1
#endif

int ctx_clean(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that we have at least 2 arguments

  // Check cmdline args
  boolean tip_cleaning = false, supernode_cleaning = false;
  size_t max_tip_len = 0;
  Covg threshold = 0;
  double seq_depth = -1;
  char *dump_covgs = NULL;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++) {
    if(strcmp(argv[argi],"--tips") == 0) {
      if(argi + 1 >= argc || !parse_entire_size(argv[argi+1], &max_tip_len) ||
         max_tip_len <= 1) {
        cmd_print_usage("--tips <L> needs an integer argument > 1");
      }
      tip_cleaning = true;
      argi++;
    }
    else if(strcmp(argv[argi],"--supernodes") == 0) supernode_cleaning = true;
    else if(strcmp(argv[argi],"--covgs") == 0) {
      if(argi + 1 >= argc) cmd_print_usage("--covgs <out.csv> needs an argument");
      dump_covgs = argv[argi+1];
      argi++;
    }
    else if(strcmp(argv[argi],"--threshold") == 0) {
      if(argi+1 >= argc || !parse_entire_uint(argv[argi+1], &threshold) ||
         threshold <= 1) {
        cmd_print_usage("--threshold <T> needs an integer argument > 1");
      }
      argi++;
    }
    else if(strcmp(argv[argi],"--kdepth") == 0) {
      if(argi+1 >= argc || !parse_entire_double(argv[argi+1], &seq_depth) ||
         seq_depth <= 1) {
        cmd_print_usage("--kdepth <C> needs a positive decimal number > 1");
      }
      argi++;
    }
    else cmd_print_usage("Unknown argument: %s", argv[argi]);
  }

  char *out_ctx_path = args->output_file_set ? args->output_file : NULL;

  if(argi == argc) cmd_print_usage("Please give input graph files");

  // default behaviour
  if(!tip_cleaning && !supernode_cleaning) {
    if(out_ctx_path != NULL)
      supernode_cleaning = tip_cleaning = true; // do both
    else
      warn("No cleaning being done: you did not specify --out <out.ctx>");
  }

  boolean doing_cleaning = (supernode_cleaning || tip_cleaning);

  if(doing_cleaning && out_ctx_path == NULL) {
    cmd_print_usage("Please specify --out <out.ctx> for cleaned graph");
  }

  if(!supernode_cleaning && threshold > 0)
    cmd_print_usage("--threshold <T> not needed if not cleaning with --supernodes");
  if(!supernode_cleaning && seq_depth > 0)
    cmd_print_usage("--kdepth <C> not needed if not cleaning with --supernodes");

  // Default behaviour
  if(supernode_cleaning && threshold != 0 && seq_depth > 0) {
    cmd_print_usage("supernode cleaning requires only one of "
                             "--threshold <T>, --depth <D>");
  }

  // Use remaining args as graph files
  char **paths = argv + argi;
  size_t i, j, num_files = (size_t)(argc - argi), total_cols = 0;

  // Open graph files
  GraphFileReader files[num_files];
  uint64_t max_ctx_kmers = 0;

  for(i = 0; i < num_files; i++)
  {
    files[i] = INIT_GRAPH_READER;
    graph_file_open(&files[i], paths[i], true);

    if(files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, files[i].hdr.kmer_size);
    }

    size_t offset = total_cols;
    total_cols += graph_file_usedcols(&files[i]);
    file_filter_update_intocol(&files[i].fltr, files[i].fltr.intocol + offset);

    max_ctx_kmers = MAX2(max_ctx_kmers, files[i].hdr.num_of_kmers);
  }

  size_t use_ncols = args->use_ncols, kmer_size = files[0].hdr.kmer_size;

  // Flatten if we don't have to remember colours / output a graph
  if(!doing_cleaning)
  {
    total_cols = use_ncols = 1;
    for(i = 0; i < num_files; i++)
      file_filter_update_intocol(&files[i].fltr, 0);
  }

  if(total_cols < use_ncols) {
    warn("I only need %zu colour%s ('--ncols %zu' ignored)",
         total_cols, (total_cols != 1 ? "s" : ""), use_ncols);
    use_ncols = total_cols;
  }

  // If no arguments given we default to removing tips <= 2*kmer_size - 1
  if(tip_cleaning && max_tip_len == 0)
    max_tip_len = 2 * kmer_size - 1;

  // Warn if any files already cleaned
  size_t fromcol, intocol;
  ErrorCleaning *cleaning;

  for(i = 0; i < num_files; i++) {
    for(j = 0; j < files[i].fltr.ncols; j++) {
      fromcol = graph_file_fromcol(&files[i], j);
      cleaning = &files[i].hdr.ginfo[fromcol].cleaning;
      if(cleaning->cleaned_snodes && supernode_cleaning) {
        warn("%s:%zu already has supernode cleaning with threshold: <%zu",
             files[i].fltr.file_path.buff, fromcol,
             (size_t)cleaning->clean_snodes_thresh);
      }
      if(cleaning->cleaned_tips && tip_cleaning) {
        warn("%s:%zu already has had tip cleaned",
             files[i].fltr.file_path.buff, fromcol);
      }
    }
  }

  // Print steps
  size_t step = 0;
  status("Actions:\n");
  if(tip_cleaning)
    status("%zu. Cleaning tips shorter than %zu nodes", step++, max_tip_len);
  if(dump_covgs != NULL)
    status("%zu. Saving coverage distribution to: %s", step++, dump_covgs);
  if(supernode_cleaning && threshold > 0)
    status("%zu. Cleaning supernodes with threshold < %u", step++, threshold);
  if(supernode_cleaning && threshold == 0)
    status("%zu. Cleaning supernodes with auto-detected threshold", step++);

  //
  // Pick hash table size
  //
  boolean all_colours_loaded = (total_cols <= use_ncols);

  size_t kmers_in_hash, extra_bits_per_kmer, graph_mem;
  extra_bits_per_kmer = (sizeof(Covg) + sizeof(Edges)) * 8 * use_ncols +
                        (!all_colours_loaded) * sizeof(Edges) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        max_ctx_kmers, true, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  // Check output files are writable
  if(out_ctx_path != NULL && !futil_is_file_writable(out_ctx_path))
    cmd_print_usage("Cannot write to output: %s", out_ctx_path);

  if(dump_covgs && !futil_is_file_writable(dump_covgs))
    cmd_print_usage("Cannot write coverage distribution to: %s", dump_covgs);

  // Create db_graph
  // Load as many colours as possible
  // Use an extra set of edge to take intersections
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, use_ncols, use_ncols, kmers_in_hash);
  Edges *edge_store = calloc2(db_graph.ht.capacity * (use_ncols+!all_colours_loaded),
                              sizeof(Edges));
  db_graph.col_edges = edge_store;
  db_graph.col_covgs = calloc2(db_graph.ht.capacity * use_ncols, sizeof(Covg));

  // Load graph into a single colour
  LoadingStats stats;
  loading_stats_init(&stats);

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = false};

  // Construct cleaned graph header
  GraphFileHeader outhdr = {.version = CTX_GRAPH_FILEFORMAT,
                            .kmer_size = (uint32_t)db_graph.kmer_size,
                            .num_of_bitfields = NUM_BKMER_WORDS,
                            .num_of_cols = (uint32_t)total_cols,
                            .num_of_kmers = db_graph.ht.num_kmers,
                            .capacity = 0};

  graph_header_alloc(&outhdr, total_cols);

  // Merge info into header
  size_t gcol = 0;
  for(i = 0; i < num_files; i++) {
    for(j = 0; j < files[i].fltr.ncols; j++, gcol++) {
      fromcol = graph_file_fromcol(&files[i], j);
      intocol = graph_file_intocol(&files[i], j);
      graph_info_merge(&outhdr.ginfo[intocol], &files[i].hdr.ginfo[fromcol]);
    }
  }

  if(total_cols > use_ncols)
  {
    // Load into one colour
    size_t tmpinto; boolean tmpflatten;
    for(i = 0; i < num_files; i++)
    {
      tmpinto = files[i].fltr.intocol; tmpflatten = files[i].fltr.flatten;
      file_filter_update_intocol(&files[i].fltr, 0);
      files[i].fltr.flatten = true;
      graph_load(&files[i], gprefs, &stats);
      file_filter_update_intocol(&files[i].fltr, tmpinto);
      files[i].fltr.flatten = tmpflatten;
    }
  }
  else {
    for(i = 0; i < num_files; i++)
      graph_load(&files[i], gprefs, &stats);
  }

  if(num_files > 1) {
    char num_kmers_str[100];
    ulong_to_str(db_graph.ht.num_kmers, num_kmers_str);
    status("Total kmers loaded: %s\n", num_kmers_str);
  }

  size_t initial_nkmers = db_graph.ht.num_kmers;

  hash_table_print_stats(&db_graph.ht);

  size_t visited_words = roundup_bits2words64(db_graph.ht.capacity);
  uint64_t *visited = calloc2(visited_words, sizeof(uint64_t));

  // Tip clipping
  if(tip_cleaning) {
    cleaning_remove_tips(max_tip_len, visited, &db_graph);

    if(supernode_cleaning || dump_covgs)
      memset(visited, 0, visited_words * sizeof(uint64_t));
  }

  // Supernode cleaning or printing coverage distribution to a file
  if(supernode_cleaning || dump_covgs) {
    threshold = cleaning_remove_supernodes(supernode_cleaning, threshold,
                                           seq_depth, dump_covgs,
                                           visited, &db_graph);
  }

  free(visited);

  if(doing_cleaning)
  {
    // Output graph file
    Edges *intersect_edges = NULL;
    boolean kmers_loaded = true;
    size_t thresh;

    // Set output header ginfo cleaned
    for(i = 0; i < total_cols; i++)
    {
      cleaning = &outhdr.ginfo[i].cleaning;
      cleaning->cleaned_snodes |= supernode_cleaning;
      cleaning->cleaned_tips |= tip_cleaning;

      if(supernode_cleaning) {
        thresh = cleaning->clean_snodes_thresh;
        thresh = cleaning->cleaned_snodes ? MIN2(thresh, threshold) : threshold;
        cleaning->clean_snodes_thresh = thresh;
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

    graph_files_merge(out_ctx_path, files, num_files,
                      kmers_loaded, all_colours_loaded,
                      intersect_edges, &outhdr, &db_graph);
  }

  ctx_check(db_graph.ht.num_kmers == hash_table_count_kmers(&db_graph.ht));

  graph_header_dealloc(&outhdr);

  free(edge_store);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);

  return EXIT_SUCCESS;
}
