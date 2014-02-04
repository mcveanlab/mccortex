#include "global.h"
#include "all_tests.h"
#include "build_graph.h"
#include "correct_alignment.h"
#include "db_alignment.h"

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
  // Path data
  graph.kmer_paths = malloc2(graph.ht.capacity * sizeof(PathIndex));
  graph.path_kmer_locks = calloc2(roundup_bits2bytes(graph.ht.capacity), 1);

  memset(graph.kmer_paths, 0xff, graph.ht.capacity * sizeof(PathIndex));
  path_store_alloc(&graph.pdata, 1024, 0, ncols);

  // mutations:                            **
  char seq[] = "ATGCATGTTGACCAAATAAGTCACTGTGGGAGCCACGTAAAGCGTTCGCACCGATTTGTG";
  char mut[] =     "ATGTTGACCAAATAAGTCACTGTCCGAGCCACGTAAAGCGTTCGCACC";
  char res[] =     "ATGTTGACCAAATAAGTCACTGTGGGAGCCACGTAAAGCGTTCGCACC";

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
  read_t r1 = {.name = {.b = rname, .end = strlen(rname), .size = strlen(rname)},
               .seq = {.b = mut, .end = strlen(mut), .size = strlen(mut)},
               .qual = {.b = empty, .end = 0, .size = 0}};

  db_alignment_from_reads(&aln, &r1, NULL, 0, 0, 0, &graph);
  correct_alignment_init(&corrector, &aln, params);

  dBNodeBuffer *nbuf;
  char outstr[100];
  nbuf = correct_alignment_nxt(&corrector);
  assert(nbuf != NULL);
  db_nodes_to_str(nbuf->data, nbuf->len, &graph, outstr);
  assert(strcmp(outstr, res) == 0);

  // Next alignment should be NULL
  nbuf = correct_alignment_nxt(&corrector);
  assert(nbuf == NULL);

  correct_aln_worker_dealloc(&corrector);
  db_alignment_dealloc(&aln);

  free(graph.bktlocks);
  free(graph.col_edges);
  free(graph.col_covgs);
  free(graph.kmer_paths);
  free(graph.path_kmer_locks);

  path_store_dealloc(&graph.pdata);
  db_graph_dealloc(&graph);
}

void test_corrected_aln()
{
  test_status("[CorrectAln] Testing correct_alignment.c");
  test_correct_aln_no_paths();
}
