#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "graph_format.h"
#include "supernode.h"
#include "prune_nodes.h"

static const char usage[] =
"usage: "CMD" clean [options] <in.ctx> [in2.ctx ...]\n"
" Clean a cortex binary. Joins binaries first, if multiple given\n"
" Clips tips before doing supernode thresholding (when doing both [default]).\n"
" Options:\n"
"  --memory <mem>        Memory to use\n"
"  --nkmers <hash-size>  Kmers in the hash table (e.g. 1G ~ 1 billion)\n"
"  --ncols <colour>      Number of samples in memory at once (speedup)\n"
"  --tips <L>            Clip tips shorter than <L> kmers\n"
"  --supernodes          Remove low coverage supernode. Additional options:\n"
"    --kdepth <C>        kmer depth: (depth*(R-Kmersize+1)/R); R = read length\n"
"    --threshold <T>     Cleaning threshold, remove supnodes where [coverage < T]\n"
"  --covgs <out.csv>     Dump covg distribution to a CSV file\n"
// "  --out <out.ctx>       Save output file\n"
"\n"
" Default: --tips 2*k-1 --supernodes\n";

// size of coverage histogram 2^11 = 2048
#define HISTSIZE 2048

#ifdef CTXVERBOSE
  #define DEBUG_SUPERNODE 1
#endif

#define supernode_covg(covgs,len) supernode_covg_mean(covgs,len)
// #define supernode_covg(covgs,len) supernode_read_starts(covgs,len)

static inline size_t supernode_covg_mean(Covg *covgs, size_t len)
{
  assert(len > 0);
  size_t i, sum = 0;
  for(i = 0; i < len; i++) sum += covgs[i];
  return (sum+len/2) / len;
}

static inline void covg_histogram(hkey_t node, dBGraph *db_graph,
                                  dBNode **nodes, Covg **tmpcovgs, size_t *ncap,
                                  uint64_t *visited, uint64_t *covg_hist)
{
  size_t i, len, cap;
  Covg reads_arriving;

  if(!bitset_has(visited, node))
  {
    cap = *ncap;
    len = supernode_find(node, nodes, ncap, db_graph);
    if(cap < *ncap) *tmpcovgs = realloc2(*tmpcovgs, *ncap * sizeof(Covg));
    for(i = 0; i < len; i++) {
      bitset_set(visited, (*nodes)[i].key);
      (*tmpcovgs)[i] = db_graph->col_covgs[(*nodes)[i].key];
    }

    reads_arriving = supernode_covg(*tmpcovgs, len);
    reads_arriving = MIN2(reads_arriving, HISTSIZE-1);
    covg_hist[reads_arriving]++;
  }
}

static inline void supernode_clean(hkey_t node, dBGraph *db_graph,
                                   dBNode **nodes, Covg **tmpcovgs, size_t *ncap,
                                   uint64_t *visited, uint32_t covg_threshold)
{
  size_t i, len, cap;
  Covg reads_arriving;

  if(!bitset_has(visited, node))
  {
    cap = *ncap;
    len = supernode_find(node, nodes, ncap, db_graph);
    if(cap < *ncap) *tmpcovgs = realloc2(*tmpcovgs, *ncap * sizeof(Covg));

    for(i = 0; i < len; i++) {
      bitset_set(visited, (*nodes)[i].key);
      (*tmpcovgs)[i] = db_graph->col_covgs[(*nodes)[i].key];
    }

    reads_arriving = supernode_covg(*tmpcovgs, len);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_thresh: ", stdout);
      db_nodes_print(*nodes, len, db_graph, stdout);
      printf(" len: %zu covg: %u threshold: %u\n", len,
             reads_arriving, covg_threshold);
    #endif

    if(reads_arriving < covg_threshold) {
      prune_supernode(db_graph, *nodes, len);
    }
  }
}

static inline void clip_tip(hkey_t node, dBGraph *db_graph,
                            dBNode **nodes, size_t *ncap,
                            uint64_t *visited, size_t min_keep_len)
{
  size_t i, len;
  Edges first, last;
  int in, out;

  if(!bitset_has(visited, node))
  {
    len = supernode_find(node, nodes, ncap, db_graph);

    for(i = 0; i < len; i++) bitset_set(visited, (*nodes)[i].key);

    first = db_node_edges_union(db_graph, (*nodes)[0].key);
    last = db_node_edges_union(db_graph, (*nodes)[len-1].key);
    in = edges_get_indegree(first, (*nodes)[0].orient);
    out = edges_get_outdegree(last, (*nodes)[len-1].orient);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_clip: ", stdout);
      db_nodes_print(*nodes, len, db_graph, stdout);
      fprintf(stdout, " len: %zu junc: %i\n", len, in+out);
    #endif

    if(in+out <= 1 && len < min_keep_len)
    {
      // char tmpstr[MAX_KMER_SIZE+1];
      // for(i = 0; i < len; i++) {
      //   printf("pruning: %s\n",
      //          binary_kmer_to_str(db_node_bkmer(db_graph, (*nodes)[i]),
      //                             db_graph->kmer_size, tmpstr));
      // }
      prune_supernode(db_graph, *nodes, len);
    }
  }
}

static void dump_covg_histogram(const char *path, uint64_t *covg_hist)
{
  status("Writing covg distribution to: %s\n", path);
  size_t i, end;
  FILE *fh = fopen(path, "w");
  fprintf(fh, "Covg,Supernodes\n");
  fprintf(fh, "1,%zu\n", (size_t)covg_hist[1]);
  for(end = HISTSIZE-1; end > 1 && covg_hist[end] == 0; end--);
  for(i = 2; i <= end; i++) {
    if(covg_hist[i] > 0)
      fprintf(fh, "%zu,%zu\n", i, (size_t)covg_hist[i]);
  }
  fclose(fh);
}

static size_t calc_supcleaning_threshold(uint64_t *covgs, size_t len,
                                         double seq_depth,
                                         const dBGraph *db_graph)
{
  assert(len > 5);
  assert(db_graph->ht.unique_kmers > 0);
  size_t i, d1len = len-2, d2len = len-3, f1, f2;
  double *tmp = malloc2((d1len+d2len) * sizeof(double));
  double *delta1 = tmp, *delta2 = tmp + d1len;

  // Get sequencing depth from coverage
  uint64_t covg_sum = 0, capacity = db_graph->ht.capacity * db_graph->num_of_cols;
  for(i = 0; i < capacity; i++) covg_sum += db_graph->col_covgs[i];
  double seq_depth_est = (double)covg_sum / db_graph->ht.unique_kmers;

  status("Kmer depth before cleaning supernodes: %.2f\n", seq_depth_est);
  if(seq_depth == -1) seq_depth = seq_depth_est;
  else status("Using sequence depth argument: %f\n", seq_depth);

  size_t fallback_thresh = MAX2(1, (seq_depth+1)/2);

  // +1 to ensure covgs is never 0
  for(i = 0; i < d1len; i++) delta1[i] = (double)(covgs[i+1]+1) / (covgs[i+2]+1);

  d1len = i;
  d2len = d1len - 1;

  if(d1len <= 2) {
    status("(using fallback1)\n");
    free(tmp);
    return fallback_thresh;
  }

  // d2len is d1len-1
  for(i = 0; i < d2len; i++) delta2[i] = delta1[i] / delta1[i+1];

  for(f1 = 0; f1 < d1len && delta1[f1] >= 1; f1++);
  for(f2 = 0; f2 < d2len && delta2[f2] > 1; f2++);

  free(tmp);

  if(f1 < d1len && f1 < (seq_depth*0.75)) { status("(using f1)\n"); return f1; }
  if(f2 < d2len) { status("(using f2)\n"); return f2; }
  else { status("(using fallback1)\n"); return fallback_thresh; }
}

// If covg_threshold is zero, uses covg distribution to calculate
// Returns covg threshold used
static uint32_t clean_supernodes(dBGraph *db_graph, boolean clean,
                                 uint32_t covg_threshold, double seq_depth,
                                 char *dump_covgs,
                                 dBNode **nodes, size_t *ncap, uint64_t *visited)
{
  if(db_graph->ht.unique_kmers == 0) return covg_threshold;
  uint64_t *covg_hist;
  size_t threshold_est, visited_words = round_bits_to_words64(db_graph->ht.capacity);
  Covg *tmpcovgs = malloc2((*ncap) * sizeof(Covg));

  if(covg_threshold == 0 || dump_covgs != NULL)
  {
    // Get supernode coverages
    covg_hist = calloc2(HISTSIZE, sizeof(uint64_t));
    HASH_TRAVERSE(&db_graph->ht, covg_histogram,
                  db_graph, nodes, &tmpcovgs, ncap, visited, covg_hist);

    if(dump_covgs != NULL) dump_covg_histogram(dump_covgs, covg_hist);

    memset(visited, 0, visited_words * sizeof(uint64_t));

    // set threshold using histogram and genome size
    threshold_est = calc_supcleaning_threshold(covg_hist, HISTSIZE,
                                               seq_depth, db_graph);

    status("Recommended supernode cleaning threshold: < %zu\n", threshold_est);

    if(covg_threshold == 0)
      covg_threshold = threshold_est;

    free(covg_hist);
  }

  if(clean)
  {
    // Remove supernodes
    status("Cleaning supernodes...\n");

    if(covg_threshold <= 1)
      warn("Supernode cleaning failed, cleaning with threshold of <= 1");
    else {
      HASH_TRAVERSE(&db_graph->ht, supernode_clean,
                    db_graph, nodes, &tmpcovgs, ncap, visited,
                    covg_threshold);
      memset(visited, 0, visited_words * sizeof(uint64_t));
    }
  }

  free(tmpcovgs);

  return covg_threshold;
}

int ctx_clean(CmdArgs *args)
{
  cmd_accept_options(args, "mnco", usage);
  // cmd_require_options(args, "o", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  // Check cmdline args
  boolean tip_cleaning = false, supernode_cleaning = false;
  uint32_t max_tip_len = 0, threshold = 0;
  double seq_depth = -1;
  char *dump_covgs = NULL;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++) {
    if(strcmp(argv[argi],"--tips") == 0) {
      if(argi + 1 >= argc || !parse_entire_uint(argv[argi+1], &max_tip_len) ||
         max_tip_len <= 1) {
        print_usage(usage, "--tips <L> needs an integer argument > 1");
      }
      tip_cleaning = true;
      argi++;
    }
    else if(strcmp(argv[argi],"--supernodes") == 0) supernode_cleaning = true;
    else if(strcmp(argv[argi],"--covgs") == 0) {
      if(argi + 1 >= argc) print_usage(usage, "--covgs <out.csv> needs an argument");
      dump_covgs = argv[argi+1];
      argi++;
    }
    else if(strcmp(argv[argi],"--threshold") == 0) {
      if(argi+1 >= argc || !parse_entire_uint(argv[argi+1], &threshold) || threshold<=1)
        print_usage(usage, "--threshold <T> needs an integer argument > 1");
      argi++;
    }
    else if(strcmp(argv[argi],"--kdepth") == 0) {
      if(argi+1 >= argc || !parse_entire_double(argv[argi+1], &seq_depth) || seq_depth <= 1)
        print_usage(usage, "--kdepth <C> needs a positive decimal number > 1");
      argi++;
    }
    else print_usage(usage, "Unknown argument: %s", argv[argi]);
  }

  if(argc - argi < 2) print_usage(usage, "Please give output and input binaries");

  // default behaviour
  if(!tip_cleaning && !supernode_cleaning) {
    if(dump_covgs == NULL) {
      // print_usage(usage, "Need at least one of: --tips, --supernode, --covgs <out>");
      supernode_cleaning = tip_cleaning = true; // do both
    }
    else
      warn("No cleaning being done: you did not specify --tips or --supernodes");
  }

  if(!supernode_cleaning && threshold > 0)
    print_usage(usage, "--threshold <T> not needed if not cleaning with --supernodes");
  if(!supernode_cleaning && seq_depth != -1)
    print_usage(usage, "--kdepth <C> not needed if not cleaning with --supernodes");

  // Default behaviour
  if(supernode_cleaning && threshold != 0 && seq_depth != -1) {
    print_usage(usage, "supernode cleaning requires only one of "
                       "--threshold <T>, --depth <D>");
  }

  // Use remaining args as binaries
  char *out_ctx_path = argv[argi++];
  char **paths = argv + argi;
  size_t i, j, num_files = argc - argi, total_cols = 0;

  // // Probe binary files
  GraphFileReader files[num_files];
  uint64_t max_ctx_kmers = 0;

  for(i = 0; i < num_files; i++)
  {
    files[i] = INIT_GRAPH_READER;
    graph_file_open(&files[i], paths[i], true);

    if(files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, files[i].hdr.kmer_size);
    }

    size_t offset = total_cols;
    total_cols += graph_file_usedcols(&files[i]);
    file_filter_update_intocol(&files[i].fltr, files[i].fltr.intocol + offset);

    max_ctx_kmers = MAX2(max_ctx_kmers, files[i].hdr.num_of_kmers);
  }

  // If no arguments given we default to clipping tips <= 2*kmer_size - 1
  if(tip_cleaning && max_tip_len == 0)
    max_tip_len = 2 * files[0].hdr.kmer_size - 1;

  size_t use_ncols = args->use_ncols;
  if(total_cols < use_ncols) {
    warn("I only need %zu colour%s ('--ncols %zu' ignored)",
         total_cols, (total_cols != 1 ? "s" : ""), use_ncols);
    use_ncols = total_cols;
  }

  // Print steps
  uint32_t step = 0;
  status("Actions:\n");
  if(tip_cleaning)
    status("%u. Cleaning tips shorter than %u nodes", step++, max_tip_len);
  if(dump_covgs != NULL)
    status("%u. Saving coverage distribution to: %s", step++, dump_covgs);
  if(supernode_cleaning && threshold > 0)
    status("%u. Cleaning supernodes with threshold < %u", step++, threshold);
  if(supernode_cleaning && threshold == 0)
    status("%u. Cleaning supernodes with auto-detected threshold", step++);

  //
  // Pick hash table size
  //
  size_t kmers_in_hash, extra_bits_per_kmer;
  extra_bits_per_kmer = (sizeof(Covg)+sizeof(Edges))*8*use_ncols + sizeof(Edges)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        max_ctx_kmers, true);

  // Check output files are writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  if(dump_covgs && !test_file_writable(dump_covgs))
    print_usage(usage, "Cannot write coverage distribution to: %s", dump_covgs);

  // Create db_graph
  // Load all data into first colour and clean, then use use_ncols to take
  // intersection.  Use an extra set of edge to take intersections
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, 1, 1, kmers_in_hash);
  Edges *edge_store = calloc2(db_graph.ht.capacity * (use_ncols+1), sizeof(Edges));
  db_graph.col_edges = edge_store;
  db_graph.col_covgs = calloc2(db_graph.ht.capacity * use_ncols, sizeof(Covg));

  // Load binary into a single colour
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = false};

  // Construct cleaned binary header
  GraphFileHeader outhdr = {.version = CTX_GRAPH_FILEFORMAT,
                            .kmer_size = db_graph.kmer_size,
                            .num_of_bitfields = NUM_BKMER_WORDS,
                            .num_of_cols = total_cols,
                            .num_of_kmers = db_graph.ht.unique_kmers,
                            .capacity = 0};

  graph_header_alloc(&outhdr, total_cols);

  // Merge info into header
  uint32_t outcol = 0;
  for(i = 0; i < num_files; i++) {
    outcol += files[i].fltr.intocol;
    for(j = 0; j < files[i].fltr.ncols; j++, outcol++)
      graph_info_merge(outhdr.ginfo + outcol, files[i].hdr.ginfo + files[i].fltr.cols[j]);
  }

  // Load into one colour
  size_t tmpinto; boolean tmpflatten;
  for(i = 0; i < num_files; i++)
  {
    tmpinto = files[i].fltr.intocol; tmpflatten = files[i].fltr.flatten;
    // files[i].fltr.intocol = 0;
    file_filter_update_intocol(&files[i].fltr, 0);
    files[i].fltr.flatten = true;
    graph_load(&files[i], &prefs, stats);
    // files[i].fltr.intocol = tmpinto;
    file_filter_update_intocol(&files[i].fltr, tmpinto);
    files[i].fltr.flatten = tmpflatten;
  }

  if(num_files > 1) {
    char num_kmers_str[100];
    ulong_to_str(db_graph.ht.unique_kmers, num_kmers_str);
    status("Total kmers: %s\n", num_kmers_str);
  }

  size_t initial_nkmers = db_graph.ht.unique_kmers;

  hash_table_print_stats(&db_graph.ht);

  // Variables possibly used
  size_t ncap = 2048;
  dBNode *nodes = malloc2(ncap * sizeof(dBNode));
  size_t visited_words = round_bits_to_words64(db_graph.ht.capacity);
  uint64_t *visited = calloc2(visited_words, sizeof(uint64_t));

  char rem_kmers_str[100];
  ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
  status("Initial kmers: %s\n", rem_kmers_str);

  // Tip clipping
  if(tip_cleaning) {
    status("Clipping tips shorter than %u...\n", max_tip_len);
    HASH_TRAVERSE(&db_graph.ht, clip_tip,
                  &db_graph, &nodes, &ncap, visited, max_tip_len);
    ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
    status("Remaining kmers: %s\n", rem_kmers_str);
    memset(visited, 0, visited_words * sizeof(uint64_t));
  }

  // Supernode cleaning or dump covg
  if(supernode_cleaning || dump_covgs) {
    threshold = clean_supernodes(&db_graph, supernode_cleaning,
                                 threshold, seq_depth,
                                 dump_covgs, &nodes, &ncap, visited);
    ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
    status("Remaining kmers: %s\n", rem_kmers_str);
  }

  free(nodes);
  free(visited);

  db_graph.ginfo[0].cleaning.remv_low_cov_sups_thresh = threshold;
  db_graph.ginfo[0].cleaning.remv_low_cov_sups = supernode_cleaning;
  db_graph.ginfo[0].cleaning.tip_clipping = tip_cleaning;

  // Set output header ginfo cleaned
  for(i = 0; i < total_cols; i++)
  {
    ErrorCleaning *cleaning = &outhdr.ginfo[i].cleaning;
    cleaning->remv_low_cov_sups_thresh
      = MIN2(threshold, cleaning->remv_low_cov_sups_thresh);
    cleaning->remv_low_cov_sups
      = MIN2(supernode_cleaning, cleaning->remv_low_cov_sups);
    cleaning->tip_clipping |= tip_cleaning;
  }

  if(supernode_cleaning || tip_cleaning)
  {
    // Output graph file
    Edges *intersect_edges = NULL;
    boolean kmers_loaded = true, colours_loaded = (total_cols == 1);

    if(total_cols > 1)
    {
      intersect_edges = db_graph.col_edges;
      db_graph.col_edges += db_graph.ht.capacity;
    }

    // Print stats on removed kmers
    size_t removed_nkmers = initial_nkmers - db_graph.ht.unique_kmers;
    double removed_pct = (100.0 * removed_nkmers) / initial_nkmers;
    char removed_str[100], initial_str[100];
    ulong_to_str(removed_nkmers, removed_str);
    ulong_to_str(initial_nkmers, initial_str);
    status("Removed %s of %s (%.2f%%) kmers", removed_str, initial_str, removed_pct);

    db_graph_realloc(&db_graph, use_ncols, use_ncols);
    graph_files_merge(out_ctx_path, files, num_files,
                      kmers_loaded, colours_loaded,
                      intersect_edges, &outhdr, &db_graph);
  }

  assert(db_graph.ht.unique_kmers == hash_table_count_assigned_nodes(&db_graph.ht));

  graph_header_dealloc(&outhdr);
  seq_loading_stats_free(stats);

  free(edge_store);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);

  return EXIT_SUCCESS;
}
