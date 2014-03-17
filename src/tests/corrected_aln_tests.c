#include "global.h"
#include "all_tests.h"
#include "build_graph.h"
#include "correct_alignment.h"
#include "db_alignment.h"

static void check_correct_aln(char *mut, const char *ans,
                              dBAlignment *aln, read_t *r1,
                              CorrectAlnWorker *corrector,
                              CorrectAlnParam params,
                              const dBGraph *graph)
{
  r1->seq.b = mut;
  r1->seq.end = strlen(mut);
  r1->seq.size = r1->seq.end+1;

  db_alignment_from_reads(aln, r1, NULL, 0, 0, 0, graph, -1);
  correct_alignment_init(corrector, aln, params);

  dBNodeBuffer *nbuf;
  char outstr[100];
  nbuf = correct_alignment_nxt(corrector);
  TASSERT(nbuf != NULL);
  db_nodes_to_str(nbuf->data, nbuf->len, graph, outstr);
  TASSERT2(strcmp(outstr, ans) == 0, "Got: %s exp: %s", outstr, ans);

  // Next alignment should be NULL
  nbuf = correct_alignment_nxt(corrector);
  TASSERT(nbuf == NULL);
}

static void test_correct_aln_no_paths()
{
    // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  dBAlignment aln;
  CorrectAlnWorker corrector;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 1024);
  // Graph data
  graph.bktlocks = calloc2(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = calloc2(graph.ht.capacity * ncols, sizeof(Edges));
  graph.col_covgs = calloc2(graph.ht.capacity * ncols, sizeof(Covg));
  graph.node_in_cols = calloc2(roundup_bits2bytes(graph.ht.capacity) * ncols, 1);

  // Path data
  path_store_alloc(&graph.pstore, 1024, 0, graph.ht.capacity, ncols);
  graph.pstore.kmer_locks = calloc2(roundup_bits2bytes(graph.ht.capacity), 1);

  // mutations:                            **                 *
  char seq[] = "ATGCATGTTGACCAAATAAGTCACTGTGGGAGCCACGTAAAGCGTTCGCACCGATTTGTG";
  char mu0[] =     "ATGTTGACCAAATAAGTCACTGTCCGAGCCACGTAAAGCGTTCGCACC";
  char re0[] =     "ATGTTGACCAAATAAGTCACTGTGGGAGCCACGTAAAGCGTTCGCACC";
  //                                    v                       *
  char mu1[] =     "ATGTTGACCAAATAAGTCA" "TGTGGGAGCCACGTAAAGCGTTAGCACCGATTTGTG";
  char re1[] =     "ATGTTGACCAAATAAGTCAC""TGTGGGAGCCACGTAAAGCGTTCGCACCGATTTGTG";

  build_graph_from_str_mt(&graph, 0, seq, strlen(seq));

  // Start alignment
  db_alignment_alloc(&aln);
  correct_aln_worker_alloc(&corrector, &graph);

  // Set up alignment correction params
  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .ins_gap_min = 0, .ins_gap_max = 0,
                            .one_way_gap_traverse = true, .max_context = 10,
                            .gap_variance = 0.1, .gap_wiggle = 5};

  // Fake reads
  char empty[10] = "", rname[20] = "Example";
  read_t r1 = {.name = {.b = rname, .end = strlen(rname), .size = 10},
               .seq  = {.b = empty, .end = 0, .size = 1},
               .qual = {.b = empty, .end = 0, .size = 1}};

  check_correct_aln(mu0, re0, &aln, &r1, &corrector, params, &graph);
  check_correct_aln(mu1, re1, &aln, &r1, &corrector, params, &graph);

  correct_aln_worker_dealloc(&corrector);
  db_alignment_dealloc(&aln);

  free(graph.bktlocks);
  free(graph.node_in_cols);
  free(graph.col_edges);
  free(graph.col_covgs);

  path_store_dealloc(&graph.pstore);
  db_graph_dealloc(&graph);
}

void test_corrected_aln()
{
  test_status("[CorrectAln] Testing correct_alignment.c");
  test_correct_aln_no_paths();
}
