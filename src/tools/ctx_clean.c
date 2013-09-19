#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "graph_format.h"
#include "supernode.h"

static const char usage[] =
"usage: "CMD" clean [options] <out.ctx> <in.ctx> [in2.ctx ...]\n"
" Clean a cortex binary. Joins binaries first, if multiple given\n"
" Clips tips before doing supernode thresholding (when doing both [default]).\n"
" Options:\n"
"  -m <mem>           Memory to use\n"
"  -h <hash-size>     Kmers in the hash table (e.g. 1G ~ 1 billion)\n"
"  --tips <L>         Clip tips shorter than <L> kmers\n"
"  --supernodes       Remove low coverage supernode. Additional options:\n"
"    --kdepth <C>     kmer depth: (depth*(R-Kmersize+1)/R); R = read length\n"
"    --threshold <T>  Cleaning threshold, remove supnodes where [coverage < T]\n"
"  --covgs <out.csv>  Dump covg distribution to a CSV file\n"
"\n"
" Default: --tips 2*kmers_size --supernodes\n";

// size of coverage histogram 2^11 = 2048
#define HISTSIZE 2048

#ifdef DEBUG
  #define DEBUG_SUPERNODE 1
#endif

static inline size_t supernode_covg_mean(Covg *covgs, size_t len)
{
  size_t i, sum = 0;
  for(i = 0; i < len; i++) sum += covgs[i];
  return (sum+len/2) / len;
}

static inline void covg_histogram(hkey_t node, dBGraph *db_graph,
                                  hkey_t **nodes, Orientation **orients,
                                  Covg **tmpcovgs, size_t *ncap,
                                  uint64_t *visited, uint64_t *covg_hist)
{
  boolean cycle;
  size_t i, len, cap;
  Covg reads_arriving;

  if(!bitset_has(visited, node))
  {
    cap = *ncap;
    len = supernode_find(db_graph, node, nodes, orients, &cycle, ncap);
    if(cap < *ncap) *tmpcovgs = realloc2(*tmpcovgs, *ncap * sizeof(Covg));
    for(i = 0; i < len; i++) {
      bitset_set(visited, (*nodes)[i]);
      (*tmpcovgs)[i] = db_graph->col_covgs[(*nodes)[i]];
    }
    // reads_arriving = supernode_read_starts(*tmpcovgs, len);
    reads_arriving = supernode_covg_mean(*tmpcovgs, len);
    reads_arriving = MIN2(reads_arriving, HISTSIZE-1);
    covg_hist[reads_arriving]++;
  }
}

static inline void supernode_clean(hkey_t node, dBGraph *db_graph,
                                   hkey_t **nodes, Orientation **orients,
                                   Covg **tmpcovgs, size_t *ncap,
                                   uint64_t *visited, uint32_t covg_threshold)
{
  boolean cycle;
  size_t i, len, cap;
  Covg reads_arriving;

  if(!bitset_has(visited, node))
  {
    cap = *ncap;
    len = supernode_find(db_graph, node, nodes, orients, &cycle, ncap);
    if(cap < *ncap) *tmpcovgs = realloc2(*tmpcovgs, *ncap * sizeof(Covg));

    for(i = 0; i < len; i++) {
      bitset_set(visited, (*nodes)[i]);
      (*tmpcovgs)[i] = db_graph->col_covgs[(*nodes)[i]];
    }

    // reads_arriving = supernode_read_starts(*tmpcovgs, len);
    reads_arriving = supernode_covg_mean(*tmpcovgs, len);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_thresh: ", stdout);
      supernode_print(stdout, db_graph, *nodes, *orients, len);
      printf(" len: %zu covg: %u threshold: %u\n", len, reads_arriving, covg_threshold);
    #endif

    if(reads_arriving < covg_threshold) {
      db_graph_prune_supernode(db_graph, *nodes, len);
    }
  }
}

static inline void clip_tip(hkey_t node, dBGraph *db_graph,
                            hkey_t **nodes, Orientation **orients,
                            size_t *ncap, uint64_t *visited, size_t min_keep_len)
{
  boolean cycle;
  size_t i, len;
  Edges first, last;
  int in, out;

  if(!bitset_has(visited, node))
  {
    len = supernode_find(db_graph, node, nodes, orients, &cycle, ncap);

    for(i = 0; i < len; i++) bitset_set(visited, (*nodes)[i]);

    first = db_node_col_edges(db_graph,0,(*nodes)[0]);
    last = db_node_col_edges(db_graph,0,(*nodes)[len-1]);
    in = edges_get_indegree(first,(*orients)[0]);
    out = edges_get_outdegree(last,(*orients)[len-1]);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_clip: ", stdout);
      supernode_print(stdout, db_graph, *nodes, *orients, len);
      fprintf(stdout, " len: %zu junc: %i\n", len, in+out);
    #endif

    if(in+out <= 1 && !cycle && len < min_keep_len) {
      db_graph_prune_supernode(db_graph, *nodes, len);
    }
  }
}

static void dump_covg_histogram(const char *path, uint64_t *covg_hist)
{
  status("Writing covg distribution to: %s\n", path);
  size_t i, end;
  FILE *fh = fopen(path, "w");
  fprintf(fh, "Covg,Supernodes\n");
  for(end = HISTSIZE-1; end > 0 && covg_hist[end] == 0; end--);
  for(i = 0; i <= end; i++) fprintf(fh, "%zu,%zu\n", i, (size_t)covg_hist[i]);
  fclose(fh);
}

static size_t calc_supcleaning_threshold(uint64_t *covgs, size_t len,
                                         double seq_depth,
                                         const dBGraph *db_graph)
{
  assert(len > 5);
  size_t i, d1len=len-2, d2len = len-3, f1, f2;
  Covg *tmp = malloc2((d1len+d2len) * sizeof(Covg));
  Covg *delta1 = tmp, *delta2 = tmp + d1len*sizeof(Covg);

  // Get sequencing depth from coverage
  uint64_t covg_sum = 0, capacity = db_graph->ht.capacity;
  for(i = 0; i < capacity; i++) covg_sum += db_graph->col_covgs[i];
  double seq_depth_est = (double)covg_sum / db_graph->ht.unique_kmers;
  status("Kmer depth before cleaning supernodes: %.2f\n", seq_depth_est);
  if(seq_depth == -1) seq_depth = seq_depth_est;
  else status("Using sequence depth argument: %f\n", seq_depth);

  size_t fallback_thresh = MAX2(1, (seq_depth+1)/2);

  for(i = 0; i < d1len && covgs[i+2] > 0; i++) delta1[i] = covgs[i+1] / covgs[i+2];

  d1len = i;
  d2len = d1len - 1;
  if(d1len <= 2) { status("(using fallback1)\n"); return fallback_thresh; }

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
                                 hkey_t **nodes, Orientation **orients,
                                 size_t *ncap, uint64_t *visited)
{
  uint64_t *covg_hist;
  size_t threshold_est, visited_words;
  Covg *tmpcovgs = malloc2(*ncap * sizeof(Covg));

  if(covg_threshold == 0 || dump_covgs != NULL)
  {
    // Get supernode coverages
    covg_hist = calloc2(HISTSIZE, sizeof(uint64_t));
    HASH_TRAVERSE(&db_graph->ht, covg_histogram,
                  db_graph, nodes, orients, &tmpcovgs, ncap, visited, covg_hist);

    if(dump_covgs != NULL) dump_covg_histogram(dump_covgs, covg_hist);

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
      die("Supernode cleaning failed, cleaning with threshold of <= 1");

    visited_words = round_bits_to_words64(db_graph->ht.capacity);
    memset(visited, 0, visited_words * sizeof(uint64_t));
    HASH_TRAVERSE(&db_graph->ht, supernode_clean,
                  db_graph, nodes, orients, &tmpcovgs, ncap, visited,
                  covg_threshold);
  }

  free(tmpcovgs);

  return covg_threshold;
}

int ctx_clean(CmdArgs *args)
{
  cmd_accept_options(args, "mhg");
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
  char **binary_paths = argv + argi;
  uint32_t i, j, num_binaries = argc - argi;

  // Probe binary files
  boolean is_binary = false;
  GraphFileHeader gheader = {.capacity = 0};
  uint32_t kmer_size = 0, num_of_cols = 0, max_cols = 0;
  uint32_t ctx_num_cols[num_binaries], ctx_max_cols[num_binaries];
  uint64_t max_ctx_kmers = 0;

  for(i = 0; i < num_binaries; i++)
  {
    if(!graph_file_probe(binary_paths[i], &is_binary, &gheader))
      print_usage(usage, "Cannot read input binary file: %s", binary_paths[i]);
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", binary_paths[i]);

    if(i == 0)
      kmer_size = gheader.kmer_size;
    else if(kmer_size != gheader.kmer_size)
      print_usage(usage, "Kmer sizes don't match [%u vs %u]", kmer_size, gheader.kmer_size);

    ctx_num_cols[i] = gheader.num_of_cols;
    ctx_max_cols[i] = gheader.max_col;
    num_of_cols += gheader.num_of_cols;
    max_cols = MAX2(max_cols, gheader.num_of_cols);
    max_ctx_kmers = MAX2(max_ctx_kmers, gheader.num_of_kmers);

    printf("%s has %u colour%s\n", binary_paths[i], gheader.num_of_cols,
           gheader.num_of_cols == 1 ? "" : "s");
  }

  // If no arguments given we default to clipping tips <= 2*kmer_size
  if(tip_cleaning && max_tip_len == 0)
    max_tip_len = 2 * kmer_size;

  uint32_t load_colours[num_binaries][max_cols];
  for(i = 0; i < num_binaries; i++)
    graph_file_parse_colours(binary_paths[i], load_colours[i], ctx_max_cols[i]);

  // Print steps
  uint32_t step = 0;
  status("Actions:\n");
  if(tip_cleaning)
    status("%u. Cleaning tips shorter than %u nodes (%u bases)\n",
            step++, max_tip_len, max_tip_len + kmer_size - 1);
  if(dump_covgs != NULL)
    status("%u. Saving coverage distribution to: %s\n", step++, dump_covgs);
  if(supernode_cleaning && threshold > 0)
    status("%u. Cleaning supernodes with threshold < %u\n", step++, threshold);
  if(supernode_cleaning && threshold == 0)
    status("%u. Cleaning supernodes with auto-detected threshold\n", step++);

  //
  // Pick hash table size
  //
  size_t kmers_in_hash, extra_bits_per_kmer = sizeof(Covg)+sizeof(Edges);
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        max_ctx_kmers, true);

  // Check output files are writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  if(dump_covgs && !test_file_writable(dump_covgs))
    print_usage(usage, "Cannot write coverage distribution to: %s", dump_covgs);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity, sizeof(Covg));

  // Load binary into a single colour
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs
    = {.into_colour = 0, .merge_colours = true,
       .boolean_covgs = false,
       .must_exist_in_graph = false,
       .empty_colours = false,
       .db_graph = &db_graph};

  // Construct cleaned binary header
  GraphFileHeader tmpheader = {.capacity = 0};
  GraphFileHeader output_header = {.version = CTX_GRAPH_FILEFORMAT,
                                   .kmer_size = db_graph.kmer_size,
                                   .num_of_bitfields = NUM_BKMER_WORDS,
                                   .num_of_cols = num_of_cols,
                                   .num_of_kmers = db_graph.ht.unique_kmers,
                                   .capacity = 0};

  graph_header_alloc(&tmpheader, max_cols);
  graph_header_alloc(&output_header, num_of_cols);

  uint32_t output_colour = 0;
  for(i = 0; i < num_binaries; i++) {
    graph_load(binary_paths[i], &prefs, stats, &tmpheader);
    for(j = 0; j < ctx_num_cols[i]; j++, output_colour++)
      graph_info_merge(output_header.ginfo + output_colour, tmpheader.ginfo + j);
  }

  char num_kmers_str[100];
  ulong_to_str(db_graph.ht.unique_kmers, num_kmers_str);
  status("Loaded %s kmers\n", num_kmers_str);

  hash_table_print_stats(&db_graph.ht);

  // Variables possibly used
  size_t ncap = 2048;
  hkey_t *nodes = malloc2(ncap * sizeof(hkey_t));
  Orientation *orients = malloc2(ncap * sizeof(Orientation));
  size_t visited_words = round_bits_to_words64(db_graph.ht.capacity);
  uint64_t *visited = calloc2(visited_words, sizeof(uint64_t));

  char rem_kmers_str[100];
  ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
  status("Initial kmers: %s\n", rem_kmers_str);

  // Tip clipping
  if(tip_cleaning) {
    status("Clipping tips shorter than %u...\n", max_tip_len);
    HASH_TRAVERSE(&db_graph.ht, clip_tip,
                  &db_graph, &nodes, &orients, &ncap, visited, max_tip_len);
    ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
    status("Remaining kmers: %s\n", rem_kmers_str);
  }

  // Supernode cleaning or dump covg
  if(supernode_cleaning || dump_covgs) {
    threshold = clean_supernodes(&db_graph, supernode_cleaning,
                                 threshold, seq_depth,
                                 dump_covgs, &nodes, &orients, &ncap, visited);
    ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
    status("Remaining kmers: %s\n", rem_kmers_str);
  }

  free(nodes);
  free(orients);
  free(visited);

  db_graph.ginfo[0].cleaning.remv_low_cov_sups_thresh = threshold;
  db_graph.ginfo[0].cleaning.remv_low_cov_sups = supernode_cleaning;
  db_graph.ginfo[0].cleaning.tip_clipping = tip_cleaning;

  // Set output header ginfo cleaned
  for(i = 0; i < num_of_cols; i++)
  {
    output_header.ginfo[i].cleaning.remv_low_cov_sups_thresh = threshold;
    output_header.ginfo[i].cleaning.remv_low_cov_sups = supernode_cleaning;
    output_header.ginfo[i].cleaning.tip_clipping = tip_cleaning;
  }

  if(supernode_cleaning || tip_cleaning)
  {
    // Output graph file
    if(num_of_cols == 1) {
      graph_file_save(out_ctx_path, &db_graph, CTX_GRAPH_FILEFORMAT, NULL, 0, 1);
    }
    else
    {
      FILE *fh = fopen(out_ctx_path, "w");
      if(fh == NULL) die("Cannot open output ctx file: %s", out_ctx_path);

      size_t header_size = graph_write_header(fh, &output_header);

      graph_write_empty(&db_graph, fh, num_of_cols);

      // load, clean and dump graph one colour at a time
      prefs.must_exist_in_graph = true;

      output_colour = 0;
      for(i = 0; i < num_binaries; i++)
      {
        for(j = 0; j < num_of_cols; j++)
        {
          memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
          memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));
          graph_load_colour(binary_paths[i], &prefs, stats, load_colours[i][j]);
          fseek(fh, header_size, SEEK_SET);
          graph_file_write_colour(&db_graph, 0, output_colour, num_of_cols, fh);
          output_colour++;
        }
      }

      fclose(fh);

      // Message not printed yet for colour at a time approach
      graph_write_status(db_graph.ht.unique_kmers, num_of_cols,
                         out_ctx_path, CTX_GRAPH_FILEFORMAT);
    }
  }

  graph_header_dealloc(&gheader);
  graph_header_dealloc(&tmpheader);
  graph_header_dealloc(&output_header);

  seq_loading_stats_free(stats);

  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
