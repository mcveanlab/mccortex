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
                         size_t dist, size_t expt_nkmers,
                         char *seq)
{
  size_t num_mask_words = roundup_bits2words64(graph->ht.capacity);
  memset(mask, 0, num_mask_words*sizeof(uint64_t));

  subgraph_from_seq(graph, dist, false, 8*graph->ht.unique_kmers, mask, &seq, 1);
  TASSERT2(graph->ht.unique_kmers == expt_nkmers,
           "kmers %"PRIu64" dist %zu", graph->ht.unique_kmers, dist);
}

void test_subgraph()
{
  test_status("[subgraph] Testing subgraph...");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 19, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000);
  // Graph data
  graph.bktlocks = calloc2(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = calloc2(graph.ht.capacity * ncols, sizeof(Edges));
  graph.col_covgs = calloc2(graph.ht.capacity * ncols, sizeof(Covg));

  size_t num_mask_words = roundup_bits2words64(graph.ht.capacity);
  uint64_t *mask = calloc2(num_mask_words, sizeof(uint64_t));

  read_t r1, r2;
  seq_read_alloc(&r1);
  seq_read_alloc(&r2);

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
  TASSERT2(graph.ht.unique_kmers == 1000-19+1, "%"PRIu64" kmers", graph.ht.unique_kmers);

  // Pull out 10, 9, ... 0 bases around 2 kmers: GAGGTGGGTCCGCCTTGCGGt
  size_t dist;
  char seed[] = "GAGGTGGGTCCGCCTTGCGGt";
  for(dist = 10; dist != SIZE_MAX; dist--)
    run_subgraph(&graph, mask, dist, 2*dist+2, seed);

  // Expect 0 kmers with an empty seed
  char seed2[] = "GGTGGGTCCGCCTTGCGGt";
  run_subgraph(&graph, mask, 100, 0, seed2);

  // 2) Rebuild graph
  // Still expect 0 kmers with an empty seed
  add_to_graph(&graph, graphseq);
  run_subgraph(&graph, mask, 100, 0, seed2);

  // 3) Rebuild graph
  // Get the whole graph
  add_to_graph(&graph, graphseq);
  run_subgraph(&graph, mask, 600, 1000-19+1, seed);

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);
  free(graph.bktlocks);
  free(graph.col_edges);
  free(graph.col_covgs);
  db_graph_dealloc(&graph);
}
