#include "global.h"
#include "all_tests.h"

#include "breakpoint_caller.h"
#include "kmer_occur.h"

static void test_kmer_occur_filter()
{
  test_status("Testing kmer_occur.c");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  const size_t kmer_size = 11, ncols = 3;
  size_t i;

  // Create graph
  db_graph_alloc(&graph, kmer_size, ncols, 1, 2000);
  graph.bktlocks = ctx_calloc(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.node_in_cols = ctx_calloc(roundup_bits2bytes(graph.ht.capacity) * ncols, 1);

  //      xyz------->>>      Y         >  <         X
  // TTCGACCCGACAGGGCAACGTAGTCCGACAGGGCACAGCCCTGTCGGGGGGTGCA

  #define NUM_NODES 3
  #define NUM_READS 3

  const char *tmp[NUM_READS]
  = {
    "AACA",
    "TTCGACCCGACAGGGCAACGTAGTCCGACAGGGCACAGCCCTGTCGGGGGGTGCA",
    "TCTAGCATGTGTGTT"};

  read_t reads[NUM_READS];
  for(i = 0; i < NUM_READS; i++) {
    seq_read_alloc(&reads[i]);
    seq_read_set(&reads[i], tmp[i]);
  }

  KOGraph kograph = kograph_create(reads, NUM_READS, true, 1, &graph);

  TASSERT(kograph.nchroms == NUM_READS);
  TASSERT(kograph.koccurs != NULL);

  // Check CCCGACAGGGCAA starts at CCCGACAGGGC
  // x=CCCGACAGGGC, y=CCGACAGGGCA, z=CGACAGGGCAA
  // X=GCCCTGTCGGG, Y=TGCCCTGTCGG, Z=TTGCCCTGTCG
  dBNode nodes[NUM_NODES];
  for(i = 0; i < NUM_NODES; i++)
    nodes[i] = db_graph_find_str(&graph, &"CCCGACAGGGCAA"[i]);

  KOccurBuffer kobuf, kobuftmp;
  kmer_occur_buf_alloc(&kobuf, 16);
  kmer_occur_buf_alloc(&kobuftmp, 16);

  size_t n = kograph_filter_stretch(kograph, nodes, NUM_NODES, &kobuf, &kobuftmp);

  // Tests
  TASSERT(n == kobuf.len);
  TASSERT2(kobuf.len == 1, "got: %zu", kobuf.len);
  TASSERT(kobuf.data[0].orient == FORWARD);
  TASSERT(kobuf.data[0].chrom == 1);
  TASSERT2(kobuf.data[0].offset == 5, "offset: %zu", (size_t)kobuf.data[0].offset);

  kmer_occur_buf_dealloc(&kobuf);
  kmer_occur_buf_dealloc(&kobuftmp);

  for(i = 0; i < NUM_READS; i++) seq_read_dealloc(&reads[i]);
  kograph_free(kograph);

  ctx_free(graph.node_in_cols);
  ctx_free(graph.bktlocks);
  db_graph_dealloc(&graph);
}

static void test_find_breakpoint()
{
  test_status("Testing breakpoint calling... [missing]");
  // DEV
}

// bubble_caller_tests.c
void test_breakpoint_caller()
{
  test_kmer_occur_filter();
  test_find_breakpoint();
}
