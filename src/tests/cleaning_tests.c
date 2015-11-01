#include "global.h"
#include "all_tests.h"

#include "db_graph.h"
#include "build_graph.h"
#include "clean_graph.h"

#include <float.h>

void _test_pick_theshold()
{
  test_status("Testing graph cleaning thresholds...");

  uint64_t kmer_covg[50] =
  {0,850162595,491126976,257485953,123269745,56011040,26052551,13244708,7794102,5359446,
   4146083,3436803,2975971,2639644,2378544,2171244,1994462,1853408,1729215,1623824,
   1531549,1446237,1374893,1313321,1254029,1200727,1151012,1108859,1068353,1032062,
   998714,967214,934593,903374,877277,851058,825934,801614,780270,756232,
   735778,719226,699749,680650,665111,647841,628028,612052,597275,23671055};

  const size_t nitems = sizeof(kmer_covg) / sizeof(kmer_covg[0]);
  int thresh;

  thresh = cleaning_pick_kmer_threshold(kmer_covg, nitems,
                                        NULL, NULL, NULL, NULL);
  TASSERT2(thresh == 20, "thresh: %i", thresh);
}

void _test_graph_cleaning()
{
  test_status("Testing graph cleaning...");

  // Construct 1 colour graph with kmer-size=19
  dBGraph graph;
  const size_t kmer_size = 19, ncols = 1, nthreads = 2;
  size_t i;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000,
                 DBG_ALLOC_EDGES | DBG_ALLOC_COVGS | DBG_ALLOC_BKTLOCKS);

  uint8_t *visited = ctx_calloc(roundup_bits2bytes(graph.ht.capacity), 1);
  uint8_t *keep    = ctx_calloc(roundup_bits2bytes(graph.ht.capacity), 1);

  // Simple graph - 1000 bases, should all be cleaned off
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

  build_graph_from_str_mt(&graph, 0, graphseq, strlen(graphseq));
  TASSERT2(graph.ht.num_kmers == 1000-19+1,
           "%"PRIu64" kmers", graph.ht.num_kmers);

  // No change (min_tip_len must be > 1)
  clean_graph(nthreads, 0, 2, NULL, NULL, visited, keep, &graph);
  TASSERT(graph.ht.num_kmers == 1000-19+1);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));

  // No change (min_tip_len must be > 1)
  clean_graph(nthreads, 0, 1000-19+1, NULL, NULL, visited, keep, &graph);
  TASSERT(graph.ht.num_kmers == 1000-19+1);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));

  // All removed
  clean_graph(nthreads, 0, 1000-19+2, NULL, NULL, visited, keep, &graph);
  TASSERT2(graph.ht.num_kmers == 0, "%"PRIu64" kmers", graph.ht.num_kmers);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));

  // Reload first 200 bases of graph 3 times
  for(i = 0; i < 3; i++)
    build_graph_from_str_mt(&graph, 0, graphseq, 200);
  TASSERT2(graph.ht.num_kmers == 200-19+1,
           "%"PRIu64" kmers", graph.ht.num_kmers);

  // First 100 bp with two SNPs
  char tmp[] =
"GGCTACCTAACCAGATATCTCTGTATcCAGCTGCATTGTGTTTAGTCTACAACGACAGAtATCCCCTTCGACGCCCGC"
"GACCTCTCTTAACGGACGACGC";

  build_graph_from_str_mt(&graph, 0, tmp, strlen(tmp));

  size_t thresh = cleaning_get_threshold(nthreads, NULL, NULL, visited, &graph);
  clean_graph(nthreads, thresh, 0, NULL, NULL, visited, keep, &graph);
  TASSERT2(thresh > 1, "threshold: %zu", thresh);

  TASSERT2(graph.ht.num_kmers == 200-19+1, "%"PRIu64" kmers", graph.ht.num_kmers);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));

  // First 78 bp with a single SNP creating a tip 23bp -> 5kmers long
  char tmp2[] =
"GGCTACCTAACCAGATATCTCTGTATACAGCTGCATTGTGTTTAGTCTACAACGACAGAAATCCCCTTCGACGgCCGC";

  // Trim off new tip
  build_graph_from_str_mt(&graph, 0, tmp2, strlen(tmp2));
  TASSERT2(graph.ht.num_kmers == 200-19+1 + 23-19+1,
           "%"PRIu64" kmers", graph.ht.num_kmers);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));
  clean_graph(nthreads, 0, 2*19-1, NULL, NULL, visited, keep, &graph);
  TASSERT2(graph.ht.num_kmers == 200-19+1, "%"PRIu64" kmers", graph.ht.num_kmers);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));

  // clear hash table + graph
  hash_table_empty(&graph.ht);
  memset(graph.col_edges, 0, ncols*graph.ht.capacity*sizeof(Edges));
  memset(graph.col_covgs, 0, ncols*graph.ht.capacity*sizeof(Covg));

  // Build a graph with a single kmer and delete it
  char tmp3[] = "AGATGTGGTTCACGGCTAG";
  build_graph_from_str_mt(&graph, 0, tmp3, strlen(tmp3));
  TASSERT2(graph.ht.num_kmers == 1, "%zu", (size_t)graph.ht.num_kmers);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));
  clean_graph(nthreads, 0, 2*19-1, NULL, NULL, visited, keep, &graph);
  TASSERT(graph.ht.num_kmers == 0, "%"PRIu64" kmers", graph.ht.num_kmers);
  TASSERT(graph.ht.num_kmers == hash_table_count_kmers(&graph.ht));

  ctx_free(visited);
  ctx_free(keep);

  db_graph_dealloc(&graph);
}

void test_cleaning()
{
  _test_pick_theshold();
  _test_graph_cleaning();
}

