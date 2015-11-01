#include "global.h"
#include "all_tests.h"

#include "kmer_occur.h"

static void test_kmer_occur_filter()
{
  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  const size_t kmer_size = 11, ncols = 3;
  size_t i;

  // Create graph
  db_graph_alloc(&graph, kmer_size, ncols, 1, 2000,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL | DBG_ALLOC_BKTLOCKS);

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

  KOGraph kograph = kograph_create(reads, NUM_READS, true, 0, 1, &graph);

  TASSERT(kograph.nchroms == NUM_READS);
  TASSERT(kograph.koccurs != NULL);

  KOccurRunBuffer koruns, koruns_tmp, koruns_ended;
  korun_buf_alloc(&koruns, 16);
  korun_buf_alloc(&koruns_tmp, 16);
  korun_buf_alloc(&koruns_ended, 16);

  // Check CCCGACAGGGCAA starts at CCCGACAGGGC
  // x=CCCGACAGGGC, y=CCGACAGGGCA, z=CGACAGGGCAA
  // X=GCCCTGTCGGG, Y=TGCCCTGTCGG, Z=TTGCCCTGTCG
  dBNode nodes[NUM_NODES];
  for(i = 0; i < NUM_NODES; i++)
    nodes[i] = db_graph_find_str(&graph, &"CCCGACAGGGCAA"[i]);

  korun_buf_reset(&koruns);
  korun_buf_reset(&koruns_ended);
  kograph_filter_extend(&kograph, nodes, NUM_NODES, true, 0, 0,
                        &koruns, &koruns_tmp, &koruns_ended);

  // Checks
  TASSERT2(koruns.len == 1, "koruns.len: %zu", koruns.len);
  TASSERT(koruns.b[0].strand == STRAND_PLUS); // left-to-right with ref
  TASSERT2(koruns.b[0].chrom == 1, "chrom: %zu", (size_t)koruns.b[0].chrom);
  TASSERT2(koruns.b[0].first == 5, "offset: %zu", (size_t)koruns.b[0].first);
  TASSERT2(koruns.b[0].last == 7, "last: %zu", (size_t)koruns.b[0].last);

  // Test reverse
  db_nodes_reverse_complement(nodes, NUM_NODES);

  korun_buf_reset(&koruns);
  korun_buf_reset(&koruns_ended);
  kograph_filter_extend(&kograph, nodes, 1, true, 0, 0, &koruns, &koruns_tmp, &koruns_ended);
  kograph_filter_extend(&kograph, nodes+1, 1, true, 0, 1, &koruns, &koruns_tmp, &koruns_ended);
  kograph_filter_extend(&kograph, nodes+2, 1, true, 0, 2, &koruns, &koruns_tmp, &koruns_ended);

  // Print out for debugging
  // printf("koruns: ");
  // koruns_print(koruns.b, koruns.len, kmer_size, stdout);
  // printf("\nkoruns_ended: ");
  // koruns_print(koruns_ended.b, koruns_ended.len, kmer_size, stdout);
  // printf("\n");

  // Check results match:
  // koruns: chromid:1:17-5:-, chromid:1:37-47:+
  // koruns_ended: chromid:1:34-24:-
  TASSERT2(koruns.len == 2, "koruns.len: %zu", koruns.len);
  TASSERT2(koruns_ended.len == 1, "koruns_ended.len: %zu", koruns_ended.len);
  TASSERT(koruns.b[0].strand == STRAND_MINUS); // reverse complement of ref
  TASSERT2(koruns.b[0].chrom == 1, "chrom: %zu", (size_t)koruns.b[0].chrom);
  TASSERT2(koruns.b[0].first == 7, "offset: %zu", (size_t)koruns.b[0].first);
  TASSERT2(koruns.b[0].last == 5, "last: %zu", (size_t)koruns.b[0].last);

  korun_buf_dealloc(&koruns);
  korun_buf_dealloc(&koruns_tmp);
  korun_buf_dealloc(&koruns_ended);

  for(i = 0; i < NUM_READS; i++) seq_read_dealloc(&reads[i]);
  kograph_dealloc(&kograph);

  db_graph_dealloc(&graph);
}

void test_kmer_occur()
{
  test_status("Testing KOGraph...");
  test_kmer_occur_filter();
}
