#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "graph_format.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
#include "graph_walker.h"
#include "repeat_walker.h"

#define DEFAULT_NUM_REPEATS 10000
#define DEFAULT_MAX_AB_DIST 1000

const char exp_abc_usage[] =
"usage: "CMD" exp_abc [options] <in.ctx>\n"
"\n"
"  Experiment in graph traversal. It measures how of the following statement holds:\n"
"    If A->B and A->B->C then B->C\n"
"  Where X->Y means if we traverse from A, we reach B\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -p, --paths <in.ctp>    Load path file (can specify multiple times)\n"
"  -N, --repeat <N>        Sample N kmers (Default "QUOTE_MACRO(DEFAULT_NUM_REPEATS)")\n"
"  -M, --max-AB-dist <M>   Test 2: Max A->B contig (Default "QUOTE_MACRO(DEFAULT_MAX_AB_DIST)")\n"
"  -P, --print             Print failed contigs\n"
"\n";

// Investigating failed contigs with:
//  ctx31 exp_abc -M 10 -N 1000 -P -p k31/perf.k31.se.ctp.gz k31/perf.k31.ctx > tmp.fa
//  ctx31 subgraph -f --seq sample.fa -dist 2 --out sample.k31.ctx k31/perf.k31.ctx
//  ~/ninja-cortex/scripts/cortex_to_graphviz.pl -k 31 --mark sample.fa -p <(zcat k31/perf.k31.se.ctp.gz) sample.k31.ctx > sample.dot
//  dot -Tpdf sample.dot > sample.pdf

static struct option longopts[] =
{
  {"help",         no_argument,       NULL, 'h'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"paths",        required_argument, NULL, 'p'},
  {"repeat",       required_argument, NULL, 'N'},
  {"max-AB-dist",  required_argument, NULL, 'M'},
  {"print",        no_argument,       NULL, 'P'},
  {NULL, 0, NULL, 0}
};

#define NUM_RESULT_VALUES 8
#define RES_ABC_SUCCESS  0
#define RES_AB_WRONG     1 /* after B->A, A->B went wrong */      /* unused */
#define RES_AB_FAILED    2 /* after B->A, couldn't get A->B */    /* unused */
#define RES_BC_WRONG     3 /* after A->B->C, B->C went wrong */
#define RES_BC_FAILED    4 /* after A->B->C, couldn't get B->C */ /* unused */
#define RES_BC_OVERSHOT  5 /* after A->B->C, B->C went too far */
#define RES_LOST_IN_RPT  6 /* All cases where we get lost in a repeat */
#define RES_NO_TRAVERSAL 7 /* Can't get anywhere from B */

typedef struct {
  size_t colour, threadid, nthreads;
  bool prime_AB; // prime the distance A->B instead of traversing
  size_t num_tests, num_limit; // Counting how many tests we've run / limit
  size_t max_AB_dist; // Max contig to assemble finding A from B
  size_t results[NUM_RESULT_VALUES];
  size_t ab_fail_state[GRPHWLK_NUM_STATES], bc_fail_state[GRPHWLK_NUM_STATES];
  dBNodeBuffer nbuf;
  GraphWalker gwlk;
  RepeatWalker rptwlk;
  const dBGraph *db_graph;
  bool print_failed_contigs;
} ExpABCWorker;

static inline void reset(GraphWalker *wlk, RepeatWalker *rptwlk,
                         const dBNodeBuffer *nbuf)
{
  graph_walker_finish(wlk);
  rpt_walker_fast_clear(rptwlk, nbuf->data, nbuf->len);
}

#define CONFIRM_SUCCESS  0
#define CONFIRM_REPEAT   1
#define CONFIRM_OVERSHOT 2
#define CONFIRM_WRONG    3
#define CONFIRM_SHORT    4

// Check we can walk along a set of nodes through the graph
// If @allow_extend is true, traverse past the end of the buffer and add nodes
static inline int confirm_seq(size_t startidx, bool allow_extend,
                              GraphWalker *wlk, RepeatWalker *rpt,
                              dBNodeBuffer *nbuf, size_t colour,
                              const dBGraph *db_graph)
{
  ctx_assert(startidx < nbuf->len);
  size_t i, init_len = nbuf->len;
  graph_walker_init(wlk, db_graph, colour, colour, nbuf->data[startidx]);

  for(i = startidx+1; graph_walker_next(wlk); i++) {
    if(!rpt_walker_attempt_traverse(rpt, wlk)) {
      reset(wlk,rpt,nbuf);
      return CONFIRM_REPEAT;
    }
    if(i < init_len) {
      if(!db_nodes_are_equal(nbuf->data[i], wlk->node)) {
        reset(wlk,rpt,nbuf);
        return CONFIRM_WRONG;
      }
    }
    else {
      db_node_buf_add(nbuf, wlk->node);
      if(!allow_extend) {
        reset(wlk,rpt,nbuf);
        nbuf->len--; // Remove node we added
        return CONFIRM_OVERSHOT;
      }
    }
  }

  // printf("stopped %zu / %zu %zu\n", i, init_len, nbuf->len);

  reset(wlk,rpt,nbuf);
  return i < init_len ? CONFIRM_SHORT : CONFIRM_SUCCESS;
}

static void print_failed(dBNode node, const dBNodeBuffer *nbuf,
                         const dBGraph *db_graph, bool is_AB, bool prime_AB)
{
  const size_t kmer_size = db_graph->kmer_size;
  char bkmerstr[MAX_KMER_SIZE+1];
  BinaryKmer bkmer = db_node_get_bkmer(db_graph, node.key);
  binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
  printf(">%s:%i %s %s\n", bkmerstr, node.orient,
         is_AB ? "A->B" : "B->C", prime_AB ? "prime_AB" : "walk_AB");
  db_nodes_print(nbuf->data, nbuf->len, db_graph, stdout);
  fputc('\n', stdout);
}

static inline
int test_statement_node(dBNode node, ExpABCWorker *wrkr)
{
  const dBGraph *db_graph = wrkr->db_graph;
  dBNodeBuffer *nbuf = &wrkr->nbuf;
  GraphWalker *wlk = &wrkr->gwlk;
  RepeatWalker *rpt = &wrkr->rptwlk;
  size_t b_idx, col = wrkr->colour;

  // rpt_walker_clear(rpt);

  db_node_buf_reset(nbuf);
  db_node_buf_add(nbuf, node);

  // size_t AB_limit = wrkr->prime_AB ? SIZE_MAX : wrkr->max_AB_dist;
  size_t walk_limit = wrkr->max_AB_dist;
  // status("walk_limit: %zu", walk_limit);

  // Walk from B to find A
  graph_walker_init(wlk, db_graph, col, col, nbuf->data[0]);
  while(graph_walker_next(wlk) && nbuf->len < walk_limit) {
    if(!rpt_walker_attempt_traverse(rpt, wlk)) {
      reset(wlk,rpt,nbuf); return RES_LOST_IN_RPT;
    }
    db_node_buf_add(nbuf, wlk->node);
  }

  reset(wlk,rpt,nbuf);

  if(nbuf->len == 1) return RES_NO_TRAVERSAL;

  // Traverse A->B
  db_nodes_reverse_complement(nbuf->data, nbuf->len);
  b_idx = nbuf->len - 1;

  if(wrkr->prime_AB)
  {
    // Prime A->B without attempting to cross
    graph_walker_prime(wlk, nbuf->data, nbuf->len, nbuf->len, true,
                       col, col, db_graph);

    while(graph_walker_next(wlk)) {
      if(!rpt_walker_attempt_traverse(rpt, wlk)) {
        reset(wlk,rpt,nbuf); return RES_LOST_IN_RPT;
      }
      db_node_buf_add(nbuf, wlk->node);
    }
  }
  else
  {
    // Attempt to traverse A->B then extend past B
    int r = confirm_seq(0, true, wlk, rpt, nbuf, col, db_graph);
    switch(r) {
      case CONFIRM_REPEAT: return RES_LOST_IN_RPT;
      case CONFIRM_OVERSHOT: ctx_assert2(0,"Can't 'overshoot' when extending");
      case CONFIRM_WRONG: return RES_AB_WRONG;
      case CONFIRM_SHORT:
        if(wrkr->print_failed_contigs)
          print_failed(node, nbuf, db_graph, true, wrkr->prime_AB);
        wrkr->ab_fail_state[wlk->last_step.status]++;
        return RES_AB_FAILED;
    }
  }

  reset(wlk,rpt,nbuf);

  if(nbuf->len == b_idx+1) return RES_NO_TRAVERSAL; // Couldn't get past B

  // Last node is now C
  // Walk from B... record whether or not we reach C
  ctx_assert(db_nodes_are_equal(nbuf->data[b_idx], db_node_reverse(node)));

  int r = confirm_seq(b_idx, false, wlk, rpt, nbuf, col, db_graph);
  switch(r) {
    case CONFIRM_REPEAT: return RES_LOST_IN_RPT;
    case CONFIRM_OVERSHOT: return RES_BC_OVERSHOT;
    case CONFIRM_WRONG: return RES_BC_WRONG;
    case CONFIRM_SHORT:
      if(wrkr->print_failed_contigs)
        print_failed(node, nbuf, db_graph, false, wrkr->prime_AB);
      wrkr->bc_fail_state[wlk->last_step.status]++;
      return RES_BC_FAILED;
    case CONFIRM_SUCCESS: return RES_ABC_SUCCESS;
  }

  die("Shouldn't reach here: r=%i", r);
  return -1;
}

// called by run_exp_abc_thread() for each entry in hash table
static inline int test_statement_bkmer(hkey_t hkey, ExpABCWorker *wrkr)
{
  // BinaryKmer bkmer = binary_kmer_from_str("TGAGGTCAGGAGTTCAAGACCAGCCTGGGCA", 31);
  // hkey = hash_table_find(&wrkr->db_graph->ht, bkmer);
  // dBNode node = {.key = hkey, .orient = FORWARD};
  // int r = test_statement_node(node, wrkr);
  // wrkr->results[r]++;
  // wrkr->num_tests++;
  // return 1;

  int r, orient;
  for(orient = 0; orient < 2 && wrkr->num_tests < wrkr->num_limit; orient++) {
    dBNode node = (dBNode){.key = hkey, .orient = orient};
    r = test_statement_node(node, wrkr);
    wrkr->results[r]++;
    wrkr->num_tests++;
  }

  // Zero means keep iterating over hash entries, non-zero means stop
  return wrkr->num_tests < wrkr->num_limit ? 0 : 1;
}

static void run_exp_abc_thread(void *ptr)
{
  ExpABCWorker *wrkr = (ExpABCWorker*)ptr;
  const dBGraph *db_graph = wrkr->db_graph;

  // // Start from each kmer, in each direction
  HASH_ITERATE_PART(&db_graph->ht, wrkr->threadid, wrkr->nthreads,
                    test_statement_bkmer, wrkr);
}

static void run_exp_abc(const dBGraph *db_graph, bool prime_AB,
                        size_t nthreads, size_t num_repeats,
                        size_t max_AB_dist, bool print_failed_contigs)
{
  ExpABCWorker *wrkrs = ctx_calloc(nthreads, sizeof(ExpABCWorker));
  size_t i, j;

  if(max_AB_dist == 0) max_AB_dist = SIZE_MAX;

  for(i = 0; i < nthreads; i++) {
    wrkrs[i].colour = 0;
    wrkrs[i].threadid = i;
    wrkrs[i].nthreads = nthreads;
    wrkrs[i].db_graph = db_graph;
    wrkrs[i].prime_AB = prime_AB;
    wrkrs[i].num_limit = num_repeats / nthreads;
    wrkrs[i].max_AB_dist = max_AB_dist;
    wrkrs[i].print_failed_contigs = print_failed_contigs;
    db_node_buf_alloc(&wrkrs[i].nbuf, 1024);
    graph_walker_alloc(&wrkrs[i].gwlk);
    rpt_walker_alloc(&wrkrs[i].rptwlk, db_graph->ht.capacity, 22); // 4MB
  }

  util_run_threads(wrkrs, nthreads, sizeof(ExpABCWorker),
                   nthreads, run_exp_abc_thread);

  // Merge results
  size_t num_tests = 0, results[NUM_RESULT_VALUES] = {0};
  size_t ab_fail_state[GRPHWLK_NUM_STATES] = {0};
  size_t bc_fail_state[GRPHWLK_NUM_STATES] = {0};

  for(i = 0; i < nthreads; i++) {
    num_tests += wrkrs[i].num_tests;
    for(j = 0; j < NUM_RESULT_VALUES; j++) results[j] += wrkrs[i].results[j];
    for(j = 0; j < GRPHWLK_NUM_STATES; j++) ab_fail_state[j] += wrkrs[i].ab_fail_state[j];
    for(j = 0; j < GRPHWLK_NUM_STATES; j++) bc_fail_state[j] += wrkrs[i].bc_fail_state[j];
    db_node_buf_dealloc(&wrkrs[i].nbuf);
    graph_walker_dealloc(&wrkrs[i].gwlk);
    rpt_walker_dealloc(&wrkrs[i].rptwlk);
  }

  // Print results
  char nrunstr[50];
  ulong_to_str(num_tests, nrunstr);
  status("Ran %s tests with %zu threads", nrunstr, nthreads);

  const char *titles[] = {"RES_ABC_SUCCESS", "RES_AB_WRONG",
                          "RES_AB_FAILED",   "RES_BC_WRONG",
                          "RES_BC_FAILED",   "RES_BC_OVERSHOT",
                          "RES_LOST_IN_RPT", "RES_NO_TRAVERSAL"};

  util_print_nums(titles, results, NUM_RESULT_VALUES, 30);

  status("AB_FAILED:");
  graph_walker_print_state_hist(ab_fail_state);
  status("BC_FAILED:");
  graph_walker_print_state_hist(bc_fail_state);

  ctx_free(wrkrs);
}

int ctx_exp_abc(int argc, char **argv)
{
  size_t i, nthreads = 0, num_repeats = 0, max_AB_dist = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  bool print_failed_contigs = false;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 't': cmd_check(!nthreads,cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_add(&gpfiles, tmp_gpfile);
        break;
      case 'N': cmd_check(!num_repeats,cmd); num_repeats = cmd_uint32_nonzero(cmd, optarg); break;
      case 'M': cmd_check(!max_AB_dist,cmd); max_AB_dist = cmd_uint32_nonzero(cmd, optarg); break;
      case 'P': cmd_check(!print_failed_contigs,cmd); print_failed_contigs = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" bubbles -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;
  if(num_repeats == 0) num_repeats = DEFAULT_NUM_REPEATS;
  if(max_AB_dist == 0) max_AB_dist = DEFAULT_MAX_AB_DIST;

  if(print_failed_contigs && nthreads != 1) {
    warn("--print forces nthreads to be one. soz.");
    nthreads = 1;
  }

  if(optind+1 != argc) cmd_print_usage("Require exactly one input graph file (.ctx)");

  const char *ctx_path = argv[optind];

  //
  // Open Graph file
  //
  GraphFileReader gfile;
  memset(&gfile, 0, sizeof(GraphFileReader));
  graph_file_open(&gfile, ctx_path);

  size_t ncols = file_filter_into_ncols(&gfile.fltr);

  // Check only loading one colour
  if(ncols > 1) die("Only implemented for one colour currently");

  // Check graph + paths are compatible
  graphs_gpaths_compatible(&gfile, 1, gpfiles.data, gpfiles.len, -1);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  // 1 bit needed per kmer if we need to keep track of kmer usage
  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 + sizeof(GPath)*8 +
                  ncols;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        gfile.num_of_kmers, gfile.num_of_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.data, gpfiles.len, ncols, rem_mem, false);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, 1, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.data, gpfiles.len, path_mem, false, &db_graph);

  // Load the graph
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
    gpath_reader_load(&gpfiles.data[i], GPATH_DIE_MISSING_KMERS, &db_graph);
    gpath_reader_close(&gpfiles.data[i]);
  }
  gpfile_buf_dealloc(&gpfiles);

  status("");
  status("Test 1: Priming region A->B (max_AB_dist: %zu)", max_AB_dist);
  run_exp_abc(&db_graph, true, nthreads, num_repeats,
              max_AB_dist, print_failed_contigs);

  status("");
  status("Test 2: Trying to traverse A->B (max_AB_dist: %zu)", max_AB_dist);
  run_exp_abc(&db_graph, false, nthreads, num_repeats,
              max_AB_dist, print_failed_contigs);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
