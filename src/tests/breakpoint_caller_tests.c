#include "global.h"
#include "all_tests.h"

#include "breakpoint_caller.h"

// bubble_caller_tests.c
void test_breakpoint_caller()
{
  test_status("Testing breakpoint calling...");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  const size_t kmer_size = 11, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000);
  // Graph data
  graph.bktlocks = ctx_calloc(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = ctx_calloc(graph.ht.capacity * ncols, sizeof(Edges));

  //   mutations:              x
  char seq0[] = "AGGGATAAAACTCTGTACTGGATCTCCCT";
  char seq1[] = "AGGGATAAAACTCTcTACTGGATCTCCCT";

  build_graph_from_str_mt(&graph, 0, seq0, strlen(seq0));
  build_graph_from_str_mt(&graph, 0, seq1, strlen(seq1));

  // DEV: do some tests here
  test_status("Warning: missing tests for breakpoint caller.");

  ctx_free(graph.bktlocks);
  ctx_free(graph.col_edges);
  db_graph_dealloc(&graph);
}
