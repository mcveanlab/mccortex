#include "global.h"

#include "tools.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_paths.h"
#include "supernode.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "seq_reader.h"

#define DEFAULT_NCONTIGS 1000

const char contigs_usage[] =
"usage: "CMD" contigs [options] <input.ctx>\n"
"  Pull out contigs from the graph, print statistics\n"
"\n"
"  Options: [ -n <nkmers> | -p <paths.ctp> ]\n"
"   --memory <M>     How much memory to use\n"
"   --colour <c>     Pull out contigs from the given colour [default: 0]\n"
"   --ncontigs <N>   Pull out <N> contigs from random kmers [default: " QUOTE_MACRO(DEFAULT_NCONTIGS) "]\n"
"   --print          Print contigs in FASTA format\n"
"   --out <out.fa>   Write contigs to a file rather than STDOUT\n"
"   --seed <in.fa>   Use seed kmers from a file. If longer than kmer-size, only\n"
"                    use the first kmer found from each input sequence.\n"
"   --no-reseed      Do not use a seed kmer if it is used in a contig\n";

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
      bkey = bkmer_get_key(bkmer, kmer_size);
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
  size_t *lengths = malloc2(capacity * sizeof(size_t));
  size_t *junctions = malloc2(capacity * sizeof(size_t));
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
    cd->lengths = realloc2(cd->lengths, cd->capacity * sizeof(size_t));
    cd->junctions = realloc2(cd->junctions, cd->capacity * sizeof(size_t));
  }
}

static inline void contig_data_dealloc(ContigData *cd)
{
  free(cd->lengths);
  free(cd->junctions);
  db_node_buf_dealloc(&cd->nodes);
}

static void pulldown_contig(hkey_t hkey, ContigData *cd,
                            const dBGraph *db_graph, size_t colour,
                            GraphWalker *wlk, RepeatWalker *rptwlk,
                            uint64_t *visited, FILE *fout)
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
      supernode_reverse(nodes->data, nodes->len);
      hkey = nodes->data[nodes->len-1].key;
    }

    dBNode node = {.key = hkey, .orient = orient};
    graph_walker_init(wlk, db_graph, colour, colour, node);

    while(graph_traverse(wlk) && rpt_walker_attempt_traverse(rptwlk, wlk))
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
  uint64_t *const visited;
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


int ctx_contigs(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked there is at least 1 argument

  size_t i, n_rand_contigs = 0, colour = 0;
  boolean print_contigs = false, no_reseed = false;
  seq_file_t *seed_file = NULL;

  while(argc > 0 && argv[0][0] == '-') {
    if(strcmp(argv[0],"--ncontigs") == 0) {
      if(argc == 1 || !parse_entire_size(argv[1], &n_rand_contigs) ||
         n_rand_contigs == 0) {
        cmd_print_usage("--ncontigs <N> requires an integer argument [>0]");
      }
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--colour") == 0 || strcmp(argv[0],"--color") == 0) {
      if(argc == 1 || !parse_entire_size(argv[1], &colour))
        cmd_print_usage("--colour <c> requires a positive integer argument");
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--seed") == 0 || strcmp(argv[0],"--seeds") == 0) {
      if(argc == 1)
        cmd_print_usage("--seed <in.fa|fq|sam> requires an argument");
      if(seed_file != NULL) {
        cmd_print_usage("Only one --seed allowed, "
                           "try --seed <(cat file1.fa file2.fa ...)");
      }
      if((seed_file = seq_open(argv[1])) == NULL)
        die("Cannot read --seed file: %s", argv[1]);
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--no-reseed") == 0) {
      no_reseed = true;
      argv++; argc--;
    }
    else if(strcmp(argv[0],"--print") == 0) {
      print_contigs = true;
      argv++; argc--;
    }
    else cmd_print_usage("Unknown argument: %s", argv[0]);
  }

  if(argc != 1) cmd_print_usage(NULL);
  char *input_ctx_path = argv[0];

  if(seed_file != NULL && n_rand_contigs > 0)
    cmd_print_usage("Please specify one of --seed and --ncontigs");

  if(seed_file == NULL && n_rand_contigs == 0) n_rand_contigs = DEFAULT_NCONTIGS;

  //
  // Open graph file
  //
  GraphFileReader gfile = INIT_GRAPH_READER;
  graph_file_open(&gfile, input_ctx_path, true);

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
  graphs_paths_compatible(&gfile, 1, pfiles, num_pfiles);

  // Check colour specified
  if(colour >= graph_file_usedcols(&gfile)) {
    cmd_print_usage("--colour is too high (%zu > %zu)",
                 colour, graph_file_usedcols(&gfile)-1);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;
  char path_mem_str[100];

  // 1 bit needed per kmer if we need to keep track of noreseed
  bits_per_kmer = sizeof(Edges)*8 + gfile.hdr.num_of_cols + sizeof(uint64_t)*8 +
                  no_reseed;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        gfile.hdr.num_of_kmers, false, &graph_mem);

  // Paths memory
  size_t tmp_path_mem = path_files_tmp_mem_required(pfiles, num_pfiles);
  path_mem = path_max_mem + tmp_path_mem;

  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s", path_mem_str);

  // Total memory
  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  //
  // Output file if printing
  //
  FILE *fout = print_contigs ? stdout : NULL;
  if(args->output_file_set) {
    if(!print_contigs)
      warn("Ignoring --out <out> argument (maybe you forgot --print ?)");
    else if((fout = fopen(args->output_file, "w")) == NULL)
      die("Cannot open output file: %s", args->output_file);
  }

  // Allocate
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, gfile.hdr.num_of_cols, 1, kmers_in_hash);

  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);
  size_t nword64 = roundup_bits2words64(db_graph.ht.capacity);

  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.node_in_cols = calloc2(bytes_per_col*gfile.hdr.num_of_cols, sizeof(uint8_t));
  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(PathIndex));
  memset(db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(PathIndex));

  path_store_alloc(&db_graph.pdata, path_max_mem, tmp_path_mem, path_max_usedcols);

  uint64_t *visited = no_reseed ? calloc2(nword64, sizeof(uint64_t)) : NULL;

  GraphWalker wlk;
  graph_walker_alloc(&wlk);

  RepeatWalker rptwlk;
  rpt_walker_alloc(&rptwlk, db_graph.ht.capacity, 22); // 4MB

  // Load graph
  LoadingStats stats;
  loading_stats_init(&stats);

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  graph_load(&gfile, gprefs, &stats);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  paths_format_merge(pfiles, num_pfiles, false, &db_graph);

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

    read_t r1, r2;
    seq_read_alloc(&r1);
    seq_read_alloc(&r2);
    seq_parse_se_sf(seed_file, 0, &r1, &r2, parse_seed_reads, &ps);
    seq_read_dealloc(&r1);
    seq_read_dealloc(&r2);
    // seq_parse_se_sf() closes seed_file
  }

  if(args->output_file_set && print_contigs && fout != stdout)
    fclose(fout);

  status("\n");

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

  num_to_str(len_mean, 1, len_mean_str);
  num_to_str(jnc_mean, 1, jnc_mean_str);
  num_to_str(len_median, 1, len_median_str);
  num_to_str(jnc_median, 1, jnc_median_str);
  num_to_str(len_n50, 1, len_n50_str);
  num_to_str(jnc_n50, 1, jnc_n50_str);
  num_to_str(cd.min_len, 1, len_min_str);
  num_to_str(cd.min_junc, 1, jnc_min_str);
  num_to_str(cd.max_len, 1, len_max_str);
  num_to_str(cd.max_junc, 1, jnc_max_str);
  num_to_str(cd.total_len, 1, len_total_str);
  num_to_str(cd.total_junc, 1, jnc_total_str);

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
  uint64_t nsteps = cd.total_len - cd.ncontigs, ncontigends = 2*cd.ncontigs;
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

  contig_data_dealloc(&cd);

  rpt_walker_dealloc(&rptwlk);
  graph_walker_dealloc(&wlk);

  if(visited != NULL) free(visited);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free(db_graph.kmer_paths);
  path_store_dealloc(&db_graph.pdata);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&gfile);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  return EXIT_SUCCESS;
}
