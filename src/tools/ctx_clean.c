#include "global.h"

#include "tools.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "graph_format.h"
#include "supernode.h"
#include "prune_nodes.h"

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

#include "objbuf_macro.h"
create_objbuf(covg_buf,CovgBuffer,Covg)

#define supernode_covg(covgs,len) supernode_covg_mean(covgs,len)
// #define supernode_covg(covgs,len) supernode_read_starts(covgs,len)

static inline size_t supernode_covg_mean(Covg *covgs, size_t len)
{
  assert(len > 0);
  size_t i, sum = 0;
  for(i = 0; i < len; i++) sum += covgs[i];
  return (sum+len/2) / len;
}

// Mark each node in the supernode as visited
static inline size_t fetch_supernode(hkey_t node,
                                     dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                                     uint64_t *visited, const dBGraph *db_graph)
{
  size_t i;
  hkey_t hkey;

  db_node_buf_reset(nbuf);
  covg_buf_reset(cbuf);
  supernode_find(node, nbuf, db_graph);
  covg_buf_ensure_capacity(cbuf, nbuf->len);
  cbuf->len = nbuf->len;

  for(i = 0; i < nbuf->len; i++) {
    hkey = nbuf->data[i].key;
    bitset_set(visited, hkey);
    cbuf->data[i] = db_graph->col_covgs[hkey];
  }

  size_t reads_arriving = supernode_covg(cbuf->data, cbuf->len);
  return reads_arriving;
}

static inline void covg_histogram(hkey_t node,
                                  dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                                  uint64_t *visited, uint64_t *covg_hist,
                                  const dBGraph *db_graph)
{
  size_t reads_arriving;

  if(!bitset_get(visited, node))
  {
    reads_arriving = fetch_supernode(node, nbuf, cbuf, visited, db_graph);
    reads_arriving = MIN2(reads_arriving, HISTSIZE-1);
    covg_hist[reads_arriving]++;
  }
}

static inline void supernode_clean(hkey_t node,
                                   dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                                   uint64_t *visited, Covg covg_threshold,
                                   dBGraph *db_graph)
{
  size_t reads_arriving;

  if(!bitset_get(visited, node))
  {
    reads_arriving = fetch_supernode(node, nbuf, cbuf, visited, db_graph);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_thresh: ", stdout);
      db_nodes_print(nbuf->data, nbuf->len, db_graph, stdout);
      printf(" len: %zu covg: %zu threshold: %u\n", nbuf->len,
             reads_arriving, covg_threshold);
    #endif

    if(reads_arriving < covg_threshold) {
      prune_supernode(nbuf->data, nbuf->len, db_graph);
    }
  }
}

static inline void clip_tip(hkey_t node, dBNodeBuffer *nbuf,
                            uint64_t *visited, size_t min_keep_len,
                            dBGraph *db_graph)
{
  size_t i;
  Edges first, last;
  int in, out;

  if(!bitset_get(visited, node))
  {
    db_node_buf_reset(nbuf);
    supernode_find(node, nbuf, db_graph);

    for(i = 0; i < nbuf->len; i++) bitset_set(visited, nbuf->data[i].key);

    first = db_node_get_edges_union(db_graph, nbuf->data[0].key);
    last = db_node_get_edges_union(db_graph, nbuf->data[nbuf->len-1].key);
    in = edges_get_indegree(first, nbuf->data[0].orient);
    out = edges_get_outdegree(last, nbuf->data[nbuf->len-1].orient);

    #ifdef DEBUG_SUPERNODE
      fputs("sn_clip: ", stdout);
      db_nodes_print(nbuf->data, nbuf->len, db_graph, stdout);
      fprintf(stdout, " len: %zu junc: %i\n", nbuf->len, in+out);
    #endif

    if(in+out <= 1 && nbuf->len < min_keep_len)
      prune_supernode(nbuf->data, nbuf->len, db_graph);
  }
}

static void dump_covg_histogram(const char *path, uint64_t *covg_hist)
{
  status("Writing covg distribution to: %s\n", path);
  size_t i, end;
  FILE *fout = fopen(path, "w");
  fprintf(fout, "Covg,Supernodes\n");
  fprintf(fout, "1,%zu\n", (size_t)covg_hist[1]);
  for(end = HISTSIZE-1; end > 1 && covg_hist[end] == 0; end--);
  for(i = 2; i <= end; i++) {
    if(covg_hist[i] > 0)
      fprintf(fout, "%zu,%zu\n", i, (size_t)covg_hist[i]);
  }
  fclose(fout);
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
  if(seq_depth <= 0) seq_depth = seq_depth_est;
  else status("Using sequence depth argument: %f\n", seq_depth);

  size_t fallback_thresh = (size_t)MAX2(1, (seq_depth+1)/2);

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
static Covg clean_supernodes(boolean clean, Covg covg_threshold, double seq_depth,
                             char *dump_covgs,
                             dBNodeBuffer *nbuf, CovgBuffer *cbuf,
                             uint64_t *visited, dBGraph *db_graph)
{
  if(db_graph->ht.unique_kmers == 0) return covg_threshold;
  uint64_t *covg_hist;
  size_t threshold_est, visited_words = roundup_bits2words64(db_graph->ht.capacity);

  size_t covgcap = 1024;
  Covg *tmpcovgs = malloc2(covgcap * sizeof(Covg));

  if(covg_threshold == 0 || dump_covgs != NULL)
  {
    // Get supernode coverages
    covg_hist = calloc2(HISTSIZE, sizeof(uint64_t));

    HASH_ITERATE(&db_graph->ht, covg_histogram,
                 nbuf, cbuf, visited, covg_hist, db_graph);

    if(dump_covgs != NULL) dump_covg_histogram(dump_covgs, covg_hist);

    memset(visited, 0, visited_words * sizeof(uint64_t));

    // set threshold using histogram and genome size
    threshold_est = calc_supcleaning_threshold(covg_hist, HISTSIZE,
                                               seq_depth, db_graph);

    status("Recommended supernode cleaning threshold: < %zu\n", threshold_est);

    if(covg_threshold == 0)
      covg_threshold = (Covg)threshold_est;

    free(covg_hist);
  }

  if(clean)
  {
    // Remove supernodes
    status("Cleaning supernodes...\n");

    if(covg_threshold <= 1)
      warn("Supernode cleaning failed, cleaning with threshold of <= 1");
    else {
      HASH_ITERATE(&db_graph->ht, supernode_clean,
                   nbuf, cbuf, visited,
                   covg_threshold, db_graph);
      memset(visited, 0, visited_words * sizeof(uint64_t));
    }
  }

  free(tmpcovgs);

  return covg_threshold;
}

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

  if((tip_cleaning || supernode_cleaning) && out_ctx_path == NULL) {
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
  size_t kmers_in_hash, extra_bits_per_kmer, graph_mem;
  extra_bits_per_kmer = (sizeof(Covg)+sizeof(Edges))*8*use_ncols + sizeof(Edges)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        max_ctx_kmers, true, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  // Check output files are writable
  if(out_ctx_path != NULL && !futil_is_file_writable(out_ctx_path))
    cmd_print_usage("Cannot write to output: %s", out_ctx_path);

  if(dump_covgs && !futil_is_file_writable(dump_covgs))
    cmd_print_usage("Cannot write coverage distribution to: %s", dump_covgs);

  // Create db_graph
  // Load all data into first colour and clean, then use use_ncols to take
  // intersection.  Use an extra set of edge to take intersections
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, 1, 1, kmers_in_hash);
  Edges *edge_store = calloc2(db_graph.ht.capacity * (use_ncols+1), sizeof(Edges));
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
                            .num_of_kmers = db_graph.ht.unique_kmers,
                            .capacity = 0};

  graph_header_alloc(&outhdr, total_cols);

  // Merge info into header
  size_t outcol = 0;
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
    graph_load(&files[i], gprefs, &stats);
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
  dBNodeBuffer nbuf;
  CovgBuffer cbuf;
  db_node_buf_alloc(&nbuf, 2048);
  covg_buf_alloc(&cbuf, 2048);

  size_t visited_words = roundup_bits2words64(db_graph.ht.capacity);
  uint64_t *visited = calloc2(visited_words, sizeof(uint64_t));

  char rem_kmers_str[100];
  ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
  status("Initial kmers: %s\n", rem_kmers_str);

  // Tip clipping
  if(tip_cleaning) {
    // Need to use _SAFE hash traverse since we remove elements in clip_tip()
    status("Clipping tips shorter than %zu...\n", max_tip_len);

    HASH_ITERATE_SAFE(&db_graph.ht, clip_tip,
                      &nbuf, visited, max_tip_len, &db_graph);

    ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
    status("Remaining kmers: %s\n", rem_kmers_str);
    memset(visited, 0, visited_words * sizeof(uint64_t));
  }

  // Supernode cleaning or dump covg
  if(supernode_cleaning || dump_covgs) {
    threshold = clean_supernodes(supernode_cleaning, threshold, seq_depth,
                                 dump_covgs, &nbuf, &cbuf, visited,
                                 &db_graph);
    ulong_to_str(db_graph.ht.unique_kmers, rem_kmers_str);
    status("Remaining kmers: %s\n", rem_kmers_str);
  }

  db_node_buf_dealloc(&nbuf);
  covg_buf_dealloc(&cbuf);
  free(visited);

  if(supernode_cleaning || tip_cleaning)
  {
    // Output graph file
    Edges *intersect_edges = NULL;
    boolean kmers_loaded = true, colours_loaded = (total_cols == 1);

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

  free(edge_store);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);

  return EXIT_SUCCESS;
}
