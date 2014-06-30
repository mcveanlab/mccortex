#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "graph_format.h"
#include "supernode.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "seq_reader.h"
#include "gpath_reader.h"
#include "gpath_checks.h"

#define DEFAULT_NCONTIGS 1000

const char contigs_usage[] =
"usage: "CMD" contigs [options] <input.ctx>\n"
"\n"
"  Pull out contigs from the graph, print statistics\n"
"\n"
"  -m, --memory <mem>   Memory to use\n"
"  -n, --nkmers <N>     Number of hash table entries (e.g. 1G ~ 1 billion)\n"
// "  -t, --threads <T>    Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -o, --out <out.fa>   Print contigs in FASTA [default: don't print]\n"
"  -c, --colour <c>     Pull out contigs from the given colour [default: 0]\n"
"  -N, --ncontigs <N>   Pull out <N> contigs from random kmers [default: " QUOTE_MACRO(DEFAULT_NCONTIGS) "]\n"
"  -s, --seed <in.fa>   Use seed kmers from a file. If longer than kmer-size, only\n"
"                       use the first kmer found from each input sequence.\n"
"  -R, --no-reseed      Do not use a seed kmer if it is used in a contig\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  // {"threads",      required_argument, NULL, 't'},
  {"paths",        required_argument, NULL, 'p'},
// command specific
  {"seed",         required_argument, NULL, 's'},
  {"no-reseed",    no_argument,       NULL, 'R'},
  {"ncontigs",     required_argument, NULL, 'N'},
  {"colour",       required_argument, NULL, 'c'},
  {"color",        required_argument, NULL, 'c'},
  {NULL, 0, NULL, 0}
};

#define MAXPATH 5

typedef struct {
  size_t ncontigs, capacity;
  size_t total_len, total_junc;
  size_t contigs_outdegree[5];
  size_t paths_held[MAXPATH], paths_pickdup[MAXPATH], paths_counter[MAXPATH];
  size_t grphwlk_steps[8]; // 8 states in graph_walker.h
  size_t *lengths, *junctions;
  size_t min_len, max_len, min_junc, max_junc;
  double max_junc_density;
  dBNodeBuffer nodes;
  size_t nprint; // id of next contig printed
  size_t num_reseed_abort, num_seed_not_found;
} ContigData;


// Returns first kmer that is in the graph
// Or HASH_NOT_FOUND if no kmers are in the graph
hkey_t seq_reader_first_node(const read_t *r, uint8_t qcutoff, uint8_t hp_cutoff,
                             size_t colour, const dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  size_t contig_start, contig_end, search_start = 0;
  BinaryKmer bkmer, bkey;
  Nucleotide nuc;
  hkey_t node;
  size_t next_base;

  ctx_assert(db_graph->node_in_cols != NULL);
  ctx_assert(r != NULL);

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qcutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qcutoff, hp_cutoff, &search_start);

    const char *contig = r->seq.b + contig_start;
    size_t contig_len = contig_end - contig_start;

    bkmer = binary_kmer_from_str(contig, kmer_size);
    bkmer = binary_kmer_right_shift_one_base(bkmer);

    for(next_base = kmer_size-1; next_base < contig_len; next_base++)
    {
      nuc = dna_char_to_nuc(contig[next_base]);
      bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
      bkey = binary_kmer_get_key(bkmer, kmer_size);
      node = hash_table_find(&db_graph->ht, bkey);

      if(node != HASH_NOT_FOUND &&
         db_node_has_col(db_graph, node, colour)) return node;
    }
  }

  return HASH_NOT_FOUND;
}

static void contig_print_path_dist(const size_t *hist, size_t n,
                                   const char *name, size_t ncontigs)
{
  char nout_str[100];
  size_t i;

  timestamp();
  message(" %s: ", name);
  for(i = 0; i < n; i++) {
    message("\t%zu:%s [%zu%%]", i, ulong_to_str(hist[i], nout_str),
            (size_t)((100.0*hist[i])/(2.0*ncontigs)+0.5));
  }
  message("\n");
}

static void contig_grphwlk_state(const char *str, size_t n, size_t ncontigs)
{
  char nout_str[100];
  status("  %s: %s\t[ %2zu%% ]", str, ulong_to_str(n, nout_str),
         (size_t)(((100.0*n)/ncontigs)+0.5));
}

static inline void contig_data_alloc(ContigData *cd, size_t capacity)
{
  size_t *lengths = ctx_malloc(capacity * sizeof(size_t));
  size_t *junctions = ctx_malloc(capacity * sizeof(size_t));
  ContigData tmp = {.ncontigs = 0, .capacity = capacity,
                    .total_len = 0, .total_junc = 0,
                    .contigs_outdegree = {0}, .paths_held = {0},
                    .paths_pickdup = {0}, .paths_counter = {0},
                    .grphwlk_steps = {0},
                    .lengths = lengths, .junctions = junctions,
                    .min_len = SIZE_MAX, .max_len = 0,
                    .min_junc = SIZE_MAX, .max_junc = 0,
                    .max_junc_density = 0, .nprint = 0,
                    .num_reseed_abort = 0, .num_seed_not_found = 0};
  db_node_buf_alloc(&tmp.nodes, 1024);
  memcpy(cd, &tmp, sizeof(ContigData));
}

static inline void contig_data_ensure_capacity(ContigData *cd, size_t capacity)
{
  if(capacity > cd->capacity) {
    cd->capacity = roundup2pow(capacity);
    cd->lengths = ctx_realloc(cd->lengths, cd->capacity * sizeof(size_t));
    cd->junctions = ctx_realloc(cd->junctions, cd->capacity * sizeof(size_t));
  }
}

static inline void contig_data_dealloc(ContigData *cd)
{
  ctx_free(cd->lengths);
  ctx_free(cd->junctions);
  db_node_buf_dealloc(&cd->nodes);
}

static void pulldown_contig(hkey_t hkey, ContigData *cd,
                            const dBGraph *db_graph, size_t colour,
                            GraphWalker *wlk, RepeatWalker *rptwlk,
                            uint8_t *visited, FILE *fout)
{
  // Don't use a visited kmer as a seed node if --no-reseed passed
  if(visited != NULL && bitset_get(visited, hkey)) {
    cd->num_reseed_abort++;
    if(fout != NULL) {
      fprintf(fout, ">contig%zu\n\n", cd->nprint);
      cd->nprint++;
    }
    return;
  }

  dBNodeBuffer *nodes = &cd->nodes;
  Orientation orient;
  size_t njunc = 0, i;

  db_node_buf_reset(nodes);
  db_node_buf_safe_add(nodes, hkey, FORWARD);

  for(orient = 0; orient < 2; orient++)
  {
    if(orient == 1) {
      db_nodes_reverse_complement(nodes->data, nodes->len);
      hkey = nodes->data[nodes->len-1].key;
    }

    dBNode node = {.key = hkey, .orient = orient};
    graph_walker_init(wlk, db_graph, colour, colour, node);

    while(graph_walker_next(wlk) && rpt_walker_attempt_traverse(rptwlk, wlk))
    {
      db_node_buf_add(nodes, wlk->node);
      cd->grphwlk_steps[wlk->last_step.status]++;
    }

    // Grab some stats
    njunc += wlk->fork_count;
    cd->paths_held[MIN2(wlk->paths.len, MAXPATH-1)]++;
    cd->paths_pickdup[MIN2(wlk->new_paths.len, MAXPATH-1)]++;
    cd->paths_counter[MIN2(wlk->cntr_paths.len, MAXPATH-1)]++;

    // Get failed status
    cd->grphwlk_steps[wlk->last_step.status]++;
    // nloop += rptwlk->nbloom_entries;

    graph_walker_finish(wlk);
    rpt_walker_fast_clear(rptwlk, nodes->data, nodes->len);
  }

  if(fout != NULL) {
    fprintf(fout, ">contig%zu\n", cd->nprint);
    db_nodes_print(nodes->data, nodes->len, db_graph, fout);
    putc('\n', fout);
    cd->nprint++;
  }

  // If --no-seed set, mark visited nodes are visited
  // Don't use to seed another contig
  if(visited != NULL) {
    for(i = 0; i < nodes->len; i++)
      bitset_set(visited, nodes->data[i].key);
  }

  // Out degree
  size_t len = nodes->len;
  hkey_t firstnode = nodes->data[0].key;
  hkey_t firstorient = opposite_orientation(nodes->data[0].orient);
  hkey_t lastnode = nodes->data[len-1].key;
  hkey_t lastorient = nodes->data[len-1].orient;

  int outdegree;
  outdegree = edges_get_outdegree(db_graph->col_edges[firstnode], firstorient);
  cd->contigs_outdegree[outdegree]++;
  outdegree = edges_get_outdegree(db_graph->col_edges[lastnode], lastorient);
  cd->contigs_outdegree[outdegree]++;

  cd->lengths[cd->ncontigs] = len;
  cd->junctions[cd->ncontigs] = njunc;
  cd->total_len += len;
  cd->total_junc += njunc;

  cd->max_junc_density = MAX2(cd->max_junc_density, (double)njunc / len);
  cd->max_len = MAX2(cd->max_len, len);
  cd->max_junc = MAX2(cd->max_junc, njunc);
  cd->min_len = MIN2(cd->min_len, len);
  cd->min_junc = MIN2(cd->min_junc, njunc);

  cd->ncontigs++;
}

struct ParseSeeds {
  ContigData *const cd;
  const dBGraph *const db_graph;
  const size_t colour;
  GraphWalker *const wlk;
  RepeatWalker *const rptwlk;
  uint8_t *const visited;
  FILE *const fout;
};

// Pull down a contig from the first node that is found in the graph
static inline void parse_seed_read(const read_t *r, struct ParseSeeds *ps)
{
  const dBGraph *db_graph = ps->db_graph;

  if(r == NULL) return;
  hkey_t node = seq_reader_first_node(r, 0, 0, ps->colour, db_graph);

  if(node != HASH_NOT_FOUND) {
    pulldown_contig(node, ps->cd, db_graph, ps->colour,
                    ps->wlk, ps->rptwlk, ps->visited, ps->fout);
  }
  else
  {
    ps->cd->num_seed_not_found++;

    if(ps->fout != NULL) {
      fprintf(ps->fout, ">contig%zu\n\n", ps->cd->nprint);
      ps->cd->nprint++;
    }
  }
}

static inline void parse_seed_reads(read_t *r1, read_t *r2,
                                    uint8_t qoffset1, uint8_t qoffset2,
                                    void *ptr)
{
  (void)qoffset1; (void)qoffset2;

  struct ParseSeeds *ps = (struct ParseSeeds*)ptr;
  contig_data_ensure_capacity(ps->cd, ps->cd->ncontigs+2);
  parse_seed_read(r1, ps);
  parse_seed_read(r2, ps);
}


int ctx_contigs(int argc, char **argv)
{
  size_t num_of_threads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t i, n_rand_contigs = 0, colour = 0;
  bool no_reseed = false;
  seq_file_t *seed_file = NULL;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

  // Arg parsing
  char cmd[100], shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o':
        if(out_path != NULL) cmd_print_usage(NULL);
        out_path = optarg;
        break;
      case 't':
        if(num_of_threads) die("%s set twice", cmd);
        num_of_threads = cmd_uint32_nonzero(cmd, optarg);
        break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg, true);
        gpfile_buf_add(&gpfiles, tmp_gpfile);
        break;
      case 's': // --seed <in.fa>
        if(seed_file != NULL) cmd_print_usage(NULL);
        else if((seed_file = seq_open(optarg)) == NULL)
          die("Cannot read --seed file: %s", optarg);
        break;
      case 'R':
        if(no_reseed) cmd_print_usage("-R given twice");
        no_reseed = true;
        break;
      case 'N': n_rand_contigs = cmd_uint32_nonzero(cmd, optarg); break;
      case 'c': colour = cmd_uint32(cmd, optarg); break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" contigs -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults
  if(num_of_threads == 0) num_of_threads = DEFAULT_NTHREADS;
  if(seed_file == NULL && n_rand_contigs == 0) n_rand_contigs = DEFAULT_NCONTIGS;

  if(optind+1 != argc)
    cmd_print_usage("Too %s arguments", optind == argc ? "few" : "many");

  char *ctx_path = argv[optind];

  if(seed_file != NULL && n_rand_contigs > 0)
    cmd_print_usage("Please specify one of --seed and --ncontigs");

  //
  // Open Graph file
  //
  GraphFileReader gfile = INIT_GRAPH_READER;
  graph_file_open(&gfile, ctx_path, true); // true => errors are fatal
  size_t ncols = graph_file_usedcols(&gfile);

  // Check for compatibility between graph files and path files
  graphs_gpaths_compatible(&gfile, 1, gpfiles.data, gpfiles.len);

  // Check colour specified
  if(colour >= ncols) {
    cmd_print_usage("-c, --colour is too high (%zu > %zu)",
                 colour, graph_file_usedcols(&gfile)-1);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  // 1 bit needed per kmer if we need to keep track of no_reseed
  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 + sizeof(GPath)*8 +
                  ncols + no_reseed;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        gfile.num_of_kmers, gfile.num_of_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t min_path_mem = 0, max_path_mem = 0;
  gpath_reader_max_mem_req(gpfiles.data, gpfiles.len,
                           ncols, kmers_in_hash,
                           false, false, false,
                           &min_path_mem, &max_path_mem);

  // Maximise path memory
  path_mem = min_path_mem;
  if(graph_mem + path_mem < memargs.mem_to_use)
    path_mem = memargs.mem_to_use - graph_mem;

  // Don't request more than needed
  path_mem = MIN2(path_mem, max_path_mem);
  cmd_print_mem(path_mem, "paths");

  // Total memory
  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Output file if printing
  //
  FILE *fout = out_path ? futil_open_output(out_path) : NULL;

  // Allocate
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, ncols, 1, kmers_in_hash);

  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);

  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.node_in_cols = ctx_calloc(bytes_per_col*ncols, 1);

  // Paths
  if(gpfiles.len > 0) {
    // Create a path store that does not tracks path counts
    gpath_store_alloc(&db_graph.gpstore,
                      db_graph.num_of_cols, db_graph.ht.capacity,
                      path_mem, false, false);
  }

  uint8_t *visited = NULL;

  if(no_reseed)
    visited = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

  GraphWalker wlk;
  graph_walker_alloc(&wlk);

  RepeatWalker rptwlk;
  rpt_walker_alloc(&rptwlk, db_graph.ht.capacity, 22); // 4MB

  // Load graph
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  graph_load(&gfile, gprefs, &stats);
  graph_file_close(&gfile);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  for(i = 0; i < gpfiles.len; i++) {
    gpath_reader_load(&gpfiles.data[i], true, &db_graph);
    gpath_reader_close(&gpfiles.data[i]);
  }
  gpfile_buf_dealloc(&gpfiles);

  status("Traversing graph in colour %zu...", colour);

  hkey_t node;
  ContigData cd;

  if(seed_file == NULL)
  {
    contig_data_alloc(&cd, n_rand_contigs);
    while(cd.ncontigs < n_rand_contigs) {
      do { node = db_graph_rand_node(&db_graph); }
      while(!db_node_has_col(&db_graph, node, colour));
      pulldown_contig(node, &cd, &db_graph, colour, &wlk, &rptwlk, visited, fout);
    }
  }
  else
  {
    // Parse seed file
    struct ParseSeeds ps = {.cd = &cd, .db_graph = &db_graph, .colour = colour,
                            .wlk = &wlk, .rptwlk = &rptwlk, .visited = visited,
                            .fout = fout};

    contig_data_alloc(&cd, 1024);

    read_t r1;
    seq_read_alloc(&r1);
    seq_parse_se_sf(seed_file, 0, &r1, parse_seed_reads, &ps);
    seq_read_dealloc(&r1);
    seq_close(seed_file);
  }

  if(fout != NULL && fout != stdout) fclose(fout);

  status("\n");

  if(cd.ncontigs > 0)
  {
    qsort(cd.lengths, cd.ncontigs, sizeof(size_t), cmp_size);
    qsort(cd.junctions, cd.ncontigs, sizeof(size_t), cmp_size);

    // Calculate N50s
    size_t len_n50, jnc_n50;
    len_n50 = calc_N50(cd.lengths, cd.ncontigs, cd.total_len);
    jnc_n50 = calc_N50(cd.junctions, cd.ncontigs, cd.total_junc);

    // Calculate medians
    double len_median, jnc_median, len_mean, jnc_mean;
    len_median = MEDIAN(cd.lengths, cd.ncontigs);
    jnc_median = MEDIAN(cd.junctions, cd.ncontigs);
    len_mean = (double)cd.total_len / cd.ncontigs;
    jnc_mean = (double)cd.total_junc / cd.ncontigs;

    // Print number of contigs
    char ncontigs_str[90], nprint_str[90], reseed_str[90], noseedkmer_str[90];
    long_to_str(cd.ncontigs, ncontigs_str);
    long_to_str(cd.nprint, nprint_str);
    long_to_str(cd.num_reseed_abort, reseed_str);
    long_to_str(cd.num_seed_not_found, noseedkmer_str);
    status("pulled out %s contigs", ncontigs_str);
    status("printed out %s contigs", nprint_str);
    if(no_reseed) {
      status("no-reseed aborted %s times", reseed_str);
      status("seed kmer not found %s times", noseedkmer_str);
    }

    char len_min_str[90], len_max_str[90], len_total_str[90];
    char len_mean_str[90], len_median_str[90], len_n50_str[90];

    char jnc_min_str[90], jnc_max_str[90], jnc_total_str[90];
    char jnc_mean_str[90], jnc_median_str[90], jnc_n50_str[90];

    // Use ulong_to_str instead of num_to_str to get better accuracy
    // e.g. 966 instead of 1K
    ulong_to_str(len_mean, len_mean_str);
    ulong_to_str(jnc_mean, jnc_mean_str);
    ulong_to_str(len_median, len_median_str);
    ulong_to_str(jnc_median, jnc_median_str);
    ulong_to_str(len_n50, len_n50_str);
    ulong_to_str(jnc_n50, jnc_n50_str);
    ulong_to_str(cd.min_len, len_min_str);
    ulong_to_str(cd.min_junc, jnc_min_str);
    ulong_to_str(cd.max_len, len_max_str);
    ulong_to_str(cd.max_junc, jnc_max_str);
    ulong_to_str(cd.total_len, len_total_str);
    ulong_to_str(cd.total_junc, jnc_total_str);

    status("Lengths: mean: %s, median: %s, N50: %s, min: %s, max: %s, total: %s [kmers]",
           len_mean_str, len_median_str, len_n50_str, len_min_str, len_max_str, len_total_str);
    status("Junctions: mean: %s, median: %s, N50: %s, min: %s, max: %s, total: %s [out >1]",
           jnc_mean_str, jnc_median_str, jnc_n50_str, jnc_min_str, jnc_max_str, jnc_total_str);
    status("Max junction density: %.2f\n", cd.max_junc_density);
    // status("Contigs looping back to a kmer: %zu [%.2f%%]\n", nloop,
    //        (100.0 * nloop) / cd.ncontigs);

    timestamp();
    message(" Outdegree: ");
    char nout_str[90];

    for(i = 0; i <= 4; i++) {
      message("\t%zu:%s [%zu%%]", i, ulong_to_str(cd.contigs_outdegree[i], nout_str),
              (size_t)((100.0*cd.contigs_outdegree[i])/(2.0*cd.ncontigs)+0.5));
    }
    message("\n");

    contig_print_path_dist(cd.paths_held, MAXPATH, "Paths held", cd.ncontigs);
    contig_print_path_dist(cd.paths_pickdup, MAXPATH, "Paths pickdup", cd.ncontigs);
    contig_print_path_dist(cd.paths_counter, MAXPATH, "Paths counter", cd.ncontigs);

    const size_t *states = cd.grphwlk_steps;
    size_t nsteps = cd.total_len - cd.ncontigs, ncontigends = 2*cd.ncontigs;
    status("Traversal succeeded because:");
    contig_grphwlk_state("Go straight   ", states[GRPHWLK_FORWARD], nsteps);
    contig_grphwlk_state("Go colour     ", states[GRPHWLK_COLFWD], nsteps);
    contig_grphwlk_state("Go paths      ", states[GRPHWLK_USEPATH], nsteps);
    status("Traversal halted because:");
    contig_grphwlk_state("No coverage   ", states[GRPHWLK_NOCOVG], ncontigends);
    contig_grphwlk_state("No colour covg", states[GRPHWLK_NOCOLCOVG], ncontigends);
    contig_grphwlk_state("No paths      ", states[GRPHWLK_NOPATHS], ncontigends);
    contig_grphwlk_state("Paths split   ", states[GRPHWLK_SPLIT_PATHS], ncontigends);
    contig_grphwlk_state("Missing paths ", states[GRPHWLK_MISSING_PATHS], ncontigends);

    size_t njunc = states[GRPHWLK_USEPATH] + states[GRPHWLK_NOPATHS] +
                   states[GRPHWLK_SPLIT_PATHS] + states[GRPHWLK_MISSING_PATHS];

    status("Junctions:");
    contig_grphwlk_state("Paths resolved", states[GRPHWLK_USEPATH], njunc);
  }

  contig_data_dealloc(&cd);
  rpt_walker_dealloc(&rptwlk);
  graph_walker_dealloc(&wlk);
  ctx_free(visited);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
