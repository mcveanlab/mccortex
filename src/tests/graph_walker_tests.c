#include "global.h"
#include "all_tests.h"
#include "build_graph.h"
#include "generate_paths.h"
#include "graph_walker.h"

static void _check_junction_gaps(GraphWalker *wlk, size_t *exp_gaps, size_t n)
{
  GraphStep step;
  size_t idx = 0;

  while(graph_walker_next(wlk)) {
    step = wlk->last_step;
    if(graph_step_status_is_fork(step.status)) {
      TASSERT2(idx < n, "idx:%zu n:%zu", idx, n);
      if(idx < n) {
        TASSERT2(step.path_gap == exp_gaps[idx], "got:%zu vs exp:%zu",
                 step.path_gap, exp_gaps[idx]);
      }
      idx++;
    }
  }

  TASSERT2(idx == n, "Didn't see expected no. of forks %zu vs %zu", idx, n);
}

static void _test_graph_walker_test1()
{
  test_status("Testing GraphWalker...");

  /*

        2         11         11   <- number of kmers
       -a+      +-b-+      +c-
        C \  5 /  A  \ 18 / A     <- alleles
           +--+       +--+
          /    \     /    \
       -A+      +-B-+      +C-
        A         G         T
  */

  // mutations:   a               b                            c
  char seq0[] = "CCGATTAAAGGGTTACTATAGCACAGGAATGGTCTGGCCTGTAAGAAGTCCAGCTTC"; // abc
  char seq1[] = "CAGATTAAAGGGTTACTGTAGCACAGGAATGGTCTGGCCTGTAAGATGTCCAGCTTC"; // ABC
  char seq2[] = "CCGATTAAAGGGTTACTGTAGCACAGGAATGGTCTGGCCTGTAAGAAGTCCAGCTTC"; // aBc

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  const char *seqs[3] = {seq0, seq1, seq2};

  // Set up alignment correction params
  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .frag_len_min = 0, .frag_len_max = 0,
                            .one_way_gap_traverse = true, .use_end_check = true,
                            .max_context = 10,
                            .gap_variance = 0.1, .gap_wiggle = 5};

  // Construct graph and paths with first two sequences only
  all_tests_construct_graph(&graph, kmer_size, ncols, seqs, 2, params);

  GraphWalker wlk;
  graph_walker_alloc(&wlk, &graph);

  dBNode node = db_graph_find_str(&graph, "CAGATTAAAGG");
  graph_walker_start(&wlk, node);
  size_t exp_gap1[2] = {5, 18};

  _check_junction_gaps(&wlk, exp_gap1, sizeof(exp_gap1) / sizeof(exp_gap1[0]));
  graph_walker_finish(&wlk);

  // Add the third read which should disrupt expected path gap
  // 4 new paths, 0 new kmer paths
  all_tests_add_paths(&graph, seqs[2], params, 4, 0);

  graph_walker_start(&wlk, node);
  size_t exp_gap2[2] = {5, 5+11+18};

  _check_junction_gaps(&wlk, exp_gap2, sizeof(exp_gap2) / sizeof(exp_gap2[0]));
  graph_walker_finish(&wlk);

  // Done
  graph_walker_dealloc(&wlk);
  db_graph_dealloc(&graph);
}

void test_graph_walker()
{
  _test_graph_walker_test1();
}
