#include "global.h"
#include "all_tests.h"

#include "db_graph.h"
#include "build_graph.h"
#include "subgraph.h"

static void run_subgraph(dBGraph *graph, uint8_t *mask,
                         size_t dist, bool invert, bool grab_unitigs,
                         size_t expt_nkmers, char *seq, size_t len)
{
  memset(mask, 0, roundup_bits2bytes(graph->ht.capacity));

  size_t nthreads = 2;
  subgraph_from_seq(graph, nthreads, dist, invert, grab_unitigs,
                    8*hash_table_nkmers(&graph->ht), mask,
                    &seq, &len, 1);

  TASSERT2(hash_table_nkmers(&graph->ht) == expt_nkmers,
           "expected %zu kmers, got %zu; dist %zu invert: %s",
           expt_nkmers, (size_t)hash_table_nkmers(&graph->ht),
           dist, invert ? "yes" : "no");
}

static void simple_subgraph_test()
{
  // Construct 1 colour graph with kmer-size=19
  dBGraph graph;
  size_t kmer_size = 19, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000,
                 DBG_ALLOC_EDGES | DBG_ALLOC_COVGS | DBG_ALLOC_BKTLOCKS);

  uint8_t *mask = ctx_calloc(roundup_bits2bytes(graph.ht.capacity), 1);

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

  _tests_add_to_graph(&graph, graphseq, 0);
  TASSERT2(hash_table_nkmers(&graph.ht) == 1000-19+1,
           "%zu kmers", (size_t)hash_table_nkmers(&graph.ht));

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
  _tests_add_to_graph(&graph, graphseq, 0);
  run_subgraph(&graph, mask, 100, false, false, 0, seed2, strlen(seed2));

  // 3) Rebuild graph
  // Get the whole graph
  _tests_add_to_graph(&graph, graphseq, 0);
  run_subgraph(&graph, mask, 600, false, false, 1000-19+1, seed, strlen(seed));

  ctx_free(mask);
  db_graph_dealloc(&graph);
}

static void test_subgraph_unitigs()
{
  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000,
                  DBG_ALLOC_EDGES | DBG_ALLOC_COVGS | DBG_ALLOC_BKTLOCKS);

  uint8_t *mask = ctx_calloc(roundup_bits2bytes(graph.ht.capacity), 1);

  // unitig of 5 kmers with a kmer fork either side
  char seq0[] = "ATGGTGCCTAGAAGGTA";
  char seq1[] = "cTGGTGCCTAGAAGGTg";
  size_t i;

  _tests_add_to_graph(&graph, seq0, 0);
  _tests_add_to_graph(&graph, seq1, 0);

  for(i = 1; i <= 5; i++)
    run_subgraph(&graph, mask, 0, false, true, 5, seq0+i, kmer_size);

  // Wipe the graph
  run_subgraph(&graph, mask, 0, false, true, 0, NULL, 0);

  // Rebuild the graph
  // and check that we only get the ege kmer for each of the four cases
  _tests_add_to_graph(&graph, seq0, 0);
  _tests_add_to_graph(&graph, seq1, 0);
  run_subgraph(&graph, mask, 0, false, true, 1, seq0, kmer_size);
  _tests_add_to_graph(&graph, seq0, 0);
  _tests_add_to_graph(&graph, seq1, 0);
  run_subgraph(&graph, mask, 0, false, true, 1, seq1, kmer_size);
  _tests_add_to_graph(&graph, seq0, 0);
  _tests_add_to_graph(&graph, seq1, 0);
  run_subgraph(&graph, mask, 0, false, true, 1, seq0+6, kmer_size);
  _tests_add_to_graph(&graph, seq0, 0);
  _tests_add_to_graph(&graph, seq1, 0);
  run_subgraph(&graph, mask, 0, false, true, 1, seq1+6, kmer_size);

  ctx_free(mask);
  db_graph_dealloc(&graph);
}

void test_subgraph()
{
  test_status("Testing subgraph...");
  simple_subgraph_test();
  test_subgraph_unitigs();
}
