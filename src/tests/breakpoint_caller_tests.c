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
  graph.col_edges = ctx_calloc(graph.ht.capacity, sizeof(Edges));

  //      xyz------->>>      y         >  <         X
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

  KOccurRunBuffer koruns, koruns_ended;
  kmer_run_buf_alloc(&koruns, 16);
  kmer_run_buf_alloc(&koruns_ended, 16);

  // Check CCCGACAGGGCAA starts at CCCGACAGGGC
  // x=CCCGACAGGGC, y=CCGACAGGGCA, z=CGACAGGGCAA
  // X=GCCCTGTCGGG, Y=TGCCCTGTCGG, Z=TTGCCCTGTCG
  dBNode nodes[NUM_NODES];
  for(i = 0; i < NUM_NODES; i++)
    nodes[i] = db_graph_find_str(&graph, &"CCCGACAGGGCAA"[i]);

  kmer_run_buf_reset(&koruns);
  kmer_run_buf_reset(&koruns_ended);
  kograph_filter_extend(kograph, nodes, NUM_NODES, true, 0, 0, &koruns, &koruns_ended, true);

  // Checks
  TASSERT2(koruns.len == 1, "koruns.len: %zu", koruns.len);
  TASSERT(koruns.data[0].strand == STRAND_PLUS); // left-to-right with ref
  TASSERT2(koruns.data[0].chrom == 1, "chrom: %zu", (size_t)koruns.data[0].chrom);
  TASSERT2(koruns.data[0].first == 5, "offset: %zu", (size_t)koruns.data[0].first);
  TASSERT2(koruns.data[0].last == 7, "last: %zu", (size_t)koruns.data[0].last);

  // Test reverse
  db_nodes_reverse_complement(nodes, NUM_NODES);

  kmer_run_buf_reset(&koruns);
  kmer_run_buf_reset(&koruns_ended);
  kograph_filter_extend(kograph, nodes, 1, true, 0, 0, &koruns, &koruns_ended, true);
  kograph_filter_extend(kograph, nodes+1, 1, true, 0, 1, &koruns, &koruns_ended, true);
  kograph_filter_extend(kograph, nodes+2, 1, true, 0, 2, &koruns, &koruns_ended, true);

  // Print out for debugging
  // printf("koruns: ");
  // koruns_print(koruns.data, koruns.len, kmer_size, stdout);
  // printf("\nkoruns_ended: ");
  // koruns_print(koruns_ended.data, koruns_ended.len, kmer_size, stdout);
  // printf("\n");

  // Check results match:
  // koruns: chromid:1:17-5:-, chromid:1:37-47:+
  // koruns_ended: chromid:1:34-24:-
  TASSERT2(koruns.len == 2, "koruns.len: %zu", koruns.len);
  TASSERT2(koruns_ended.len == 1, "koruns_ended.len: %zu", koruns_ended.len);
  TASSERT(koruns.data[0].strand == STRAND_MINUS); // reverse complement of ref
  TASSERT2(koruns.data[0].chrom == 1, "chrom: %zu", (size_t)koruns.data[0].chrom);
  TASSERT2(koruns.data[0].first == 7, "offset: %zu", (size_t)koruns.data[0].first);
  TASSERT2(koruns.data[0].last == 5, "last: %zu", (size_t)koruns.data[0].last);

  kmer_run_buf_dealloc(&koruns);
  kmer_run_buf_dealloc(&koruns_ended);

  for(i = 0; i < NUM_READS; i++) seq_read_dealloc(&reads[i]);
  kograph_free(kograph);

  ctx_free(graph.col_edges);
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
