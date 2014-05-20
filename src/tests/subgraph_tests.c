#include "global.h"
#include "all_tests.h"

#include "db_graph.h"
#include "build_graph.h"
#include "subgraph.h"

static void add_to_graph(dBGraph *graph, const char *str)
{
  build_graph_from_str_mt(graph, 0, str, strlen(str));
}

static void run_subgraph(dBGraph *graph, uint64_t *mask,
                         size_t dist, bool invert, bool grab_supernodes,
                         size_t expt_nkmers, char *seq, size_t len)
{
  size_t num_mask_words = roundup_bits2words64(graph->ht.capacity);
  memset(mask, 0, num_mask_words*sizeof(uint64_t));

  subgraph_from_seq(graph, dist, invert, grab_supernodes,
                    8*graph->ht.num_kmers, mask, &seq, &len, 1);

  TASSERT2(graph->ht.num_kmers == expt_nkmers,
           "expected %zu kmers, got %zu; dist %zu",
           expt_nkmers, (size_t)graph->ht.num_kmers, dist);
}

static void simple_subgraph_test()
{
  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 19, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000);
  // Graph data
  graph.bktlocks = ctx_calloc(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = ctx_calloc(graph.ht.capacity * ncols, sizeof(Edges));
  graph.col_covgs = ctx_calloc(graph.ht.capacity * ncols, sizeof(Covg));

  size_t num_mask_words = roundup_bits2words64(graph.ht.capacity);
  uint64_t *mask = ctx_calloc(num_mask_words, sizeof(uint64_t));

  // Simple graph - 1000 bases
  char graphseq[] =
"GGCTACCTAACCAGATATCTCTGTATACAGCTGCATTGTGTTTAGTCTACAACGACAGAAATCCCCTTCGACGCCCGC"
"GACCTCTCTTAACGGACGACGCCTTCCGGTTGCGATATCGATGGATCGACAGAACAAGCCGCTTCCCTAACAACTGCG"
"CATGAAATCCAAAGTGCGCCGATGCTTGCTTGACGATTCCAAATCCCCATGTGACCTGTGAAGACGACTACCGTAAGA"
"TGTGTCACGGGTCAGTCGCTTTTACCACCTACGGAAGGTAGACGGTTATACTCAATTATTGGCACTTTAGCTGGGCAG"
"GTCAAAGGGAACAAGTCTGAAGTAGATATAACCTCAGTCCTTTATACGCACGTGACCCGCGTATAATCTTGCCGGTGC"
"GCAACGAGGGGCTTGGATAAAACAGCTTGGGACTTATACGTTCACCCACGACCCGCCTTAGCTCAACGCTCGTAACGA"
"CTGAATATGAGTAACGTACCTGAGGTGGGTCCGCCTTGCGGAGGTGGTGGTTCTTACTTCTATCCTCTTGTAGAGAAA"
"AGAATAGGTCGTCACTAACACTCTTGTGGGGACAAACGTGTATCGATTCCCAAACGTCCGTTAGTGAATATCCTACGT"
"GTTCCATTCGATCACACTGGAATATGGCCTTAGTTGGCCCATCTTAGTGCGCCAAGTGTTCGCAGTGGTCGTAGGCAA"
"CAGGCATCGGCGGTCTAGAGTTCACGCCAAGTCGGCCGTGTGAAGTTAAGCGTAAGTGCGGGACAACAAACCGAATGT"
"TCCGTGGCACACATGTTCGCTTATTATCAGGTAACCCTCATCTCCAGGGAGAACGCCTCAGCAGGCTTGCACCGCTTG"
"TAATCCCTCCTTATCAGAAGTAATCGTCGTTGCCGAGTTAGATCATGTCGGGACGTTGCCCTCAAGACGCCCAACGGA"
"AAAATTCACGATAGTGGCGCTCGGGAGGAGTACGCAACTCAGCACCCCGGTGAGTAGCTCCCTT";

  add_to_graph(&graph, graphseq);
  TASSERT2(graph.ht.num_kmers == 1000-19+1, "%"PRIu64" kmers", graph.ht.num_kmers);

  // Pull out 10, 9, ... 0 bases around 2 kmers: GAGGTGGGTCCGCCTTGCGGt
  size_t dist;
  char seed[] = "GAGGTGGGTCCGCCTTGCGGt";
  for(dist = 10; dist != SIZE_MAX; dist--)
    run_subgraph(&graph, mask, dist, false, false, 2*dist+2, seed, strlen(seed));

  // Expect 0 kmers with an empty seed
  char seed2[] = "GGTGGGTCCGCCTTGCGGt";
  run_subgraph(&graph, mask, 100, false, false, 0, seed2, strlen(seed2));

  // 2) Rebuild graph
  // Still expect 0 kmers with an empty seed
  add_to_graph(&graph, graphseq);
  run_subgraph(&graph, mask, 100, false, false, 0, seed2, strlen(seed2));

  // 3) Rebuild graph
  // Get the whole graph
  add_to_graph(&graph, graphseq);
  run_subgraph(&graph, mask, 600, false, false, 1000-19+1, seed, strlen(seed));

  ctx_free(mask);
  ctx_free(graph.bktlocks);
  ctx_free(graph.col_edges);
  ctx_free(graph.col_covgs);
  db_graph_dealloc(&graph);
}

static void test_subgraph_supernodes()
{
  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000);
  // Graph data
  graph.bktlocks = ctx_calloc(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = ctx_calloc(graph.ht.capacity * ncols, sizeof(Edges));
  graph.col_covgs = ctx_calloc(graph.ht.capacity * ncols, sizeof(Covg));

  size_t num_mask_words = roundup_bits2words64(graph.ht.capacity);
  uint64_t *mask = ctx_calloc(num_mask_words, sizeof(uint64_t));

  // Supernode of 5 kmers with a kmer fork either side
  char seq0[] = "ATGGTGCCTAGAAGGTA";
  char seq1[] = "cTGGTGCCTAGAAGGTg";
  size_t i;

  add_to_graph(&graph, seq0);
  add_to_graph(&graph, seq1);

  for(i = 1; i <= 5; i++)
    run_subgraph(&graph, mask, 0, false, true, 5, seq0+i, kmer_size);

  // Wipe the graph
  run_subgraph(&graph, mask, 0, false, true, 0, NULL, 0);

  // Rebuild the graph
  // and check that we only get the ege kmer for each of the four cases
  add_to_graph(&graph, seq0);
  add_to_graph(&graph, seq1);
  run_subgraph(&graph, mask, 0, false, true, 1, seq0, kmer_size);
  add_to_graph(&graph, seq0);
  add_to_graph(&graph, seq1);
  run_subgraph(&graph, mask, 0, false, true, 1, seq1, kmer_size);
  add_to_graph(&graph, seq0);
  add_to_graph(&graph, seq1);
  run_subgraph(&graph, mask, 0, false, true, 1, seq0+6, kmer_size);
  add_to_graph(&graph, seq0);
  add_to_graph(&graph, seq1);
  run_subgraph(&graph, mask, 0, false, true, 1, seq1+6, kmer_size);

  ctx_free(mask);
  ctx_free(graph.bktlocks);
  ctx_free(graph.col_edges);
  ctx_free(graph.col_covgs);
  db_graph_dealloc(&graph);
}

void test_subgraph()
{
  test_status("Testing subgraph...");
  simple_subgraph_test();
  test_subgraph_supernodes();
}
