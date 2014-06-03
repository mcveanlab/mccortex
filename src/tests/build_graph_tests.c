#include "global.h"
#include "all_tests.h"
#include "db_graph.h"
#include "db_node.h"
#include "build_graph.h"

#include <math.h>

static Covg kmer_get_covg(const char *kmer, const dBGraph *db_graph)
{
  BinaryKmer bkmer = binary_kmer_from_str(kmer, db_graph->kmer_size);
  dBNode node = db_graph_find(db_graph, bkmer);
  return db_node_get_covg(db_graph, node.key, 0);
}

void test_build_graph()
{
  test_status("Testing remove PCR duplicates in build_graph.c");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 19, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 1024);
  // Graph data
  graph.bktlocks = ctx_calloc(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = ctx_calloc(graph.ht.capacity * ncols, sizeof(Edges));
  graph.col_covgs = ctx_calloc(graph.ht.capacity * ncols, sizeof(Covg));

  // 1 bit for forward, 1 bit for reverse per kmer
  graph.readstrt = ctx_calloc(roundup_bits2bytes(graph.ht.capacity)*2,
                           sizeof(uint8_t));

  read_t r1, r2;
  seq_read_alloc(&r1);
  seq_read_alloc(&r2);

  LoadingStats stats = LOAD_STATS_INIT_MACRO;
  size_t total_seq = 0, contigs_loaded = 0;

  // Test loading empty reads are ok
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_FF,
                            &stats, 0, &graph);

  // Load a pair of reads
  seq_read_set(&r1, "CTACGATGTATGCTTAGCTGTTCCG");
  seq_read_set(&r2, "TAGAACGTTCCCTACACGTCCTATG");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_FF,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 1);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 1);
  total_seq += r1.seq.end + r2.seq.end;
  contigs_loaded += 2;

  // Check we filter out a duplicate FF
  seq_read_set(&r1, "CTACGATGTATGCTTAGCTAATGAT");
  seq_read_set(&r2, "TAGAACGTTCCCTACACGTTGTTTG");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_FF,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 1);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 1);

  // Check we filter out a duplicate FR
  // revcmp TAGAACGTTCCCTACACGT -> AGCTAAGCATACATCGTAG
  seq_read_set(&r1, "CTACGATGTATGCTTAGCTCCGAAG");
  seq_read_set(&r2, "AGACTAAGCTAAGCATACATCGTAG");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_FR,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 1);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 1);

  // Check we filter out a duplicate RF
  // revcmp CTACGATGTATGCTTAGCT -> ACGTGTAGGGAACGTTCTA
  seq_read_set(&r1, "AGGAGTTGTCTTCTAAGGAAACGTGTAGGGAACGTTCTA");
  seq_read_set(&r2, "TAGAACGTTCCCTACACGTTTTCCACGAGTTAATCTAAG");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_RF,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 1);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 1);

  // Check we filter out a duplicate RF
  // revcmp CTACGATGTATGCTTAGCT -> ACGTGTAGGGAACGTTCTA
  // revcmp TAGAACGTTCCCTACACGT -> AGCTAAGCATACATCGTAG
  seq_read_set(&r1, "AACCCTAAAAACGTGTAGGGAACGTTCTA");
  seq_read_set(&r2, "AATGCGTGTTAGCTAAGCATACATCGTAG");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_RR,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 1);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 1);

  // Check add a duplicate when filtering is turned off
  seq_read_set(&r1, "CTACGATGTATGCTTAGCTAATGAT");
  seq_read_set(&r2, "TAGAACGTTCCCTACACGTTGTTTG");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, false, READPAIR_FF,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 2);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 2);
  total_seq += r1.seq.end + r2.seq.end;
  contigs_loaded += 2;

  // Check SE duplicate removal with FF reads
  seq_read_set(&r1, "CTACGATGTATGCTTAGCTAGTGTGATATCCTCC");
  build_graph_from_reads_mt(&r1, NULL, 0, 9, 9, 9, true, READPAIR_FF,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 2);

  // Check SE duplicate removal with RR reads
  seq_read_set(&r1, "GCGTTACCTACTGACAGCTAAGCATACATCGTAG");
  build_graph_from_reads_mt(&r1, NULL, 0, 9, 9, 9, true, READPAIR_RR,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 2);

  // Check we don't filter out reads when kmers in opposite direction
  // revcmp CTACGATGTATGCTTAGCT -> ACGTGTAGGGAACGTTCTA
  // revcmp TAGAACGTTCCCTACACGT -> AGCTAAGCATACATCGTAG
  seq_read_set(&r1, "ACGTGTAGGGAACGTTCTA""CTTCTACCGGAGGAT");
  seq_read_set(&r2, "AGCTAAGCATACATCGTAG""TACAATGCACCCTCC");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_FF,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 3);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 3);
  total_seq += r1.seq.end + r2.seq.end;
  contigs_loaded += 2;

  // shouldn't work a second time
  // revcmp CTACGATGTATGCTTAGCT -> ACGTGTAGGGAACGTTCTA
  // revcmp TAGAACGTTCCCTACACGT -> AGCTAAGCATACATCGTAG
  seq_read_set(&r1, "ACGTGTAGGGAACGTTCTA""CTTCTACCGGAGGAT");
  seq_read_set(&r2, "AGCTAAGCATACATCGTAG""TACAATGCACCCTCC");
  build_graph_from_reads_mt(&r1, &r2, 0, 9, 9, 9, true, READPAIR_FF,
                            &stats, 0, &graph);
  TASSERT(kmer_get_covg("CTACGATGTATGCTTAGCT", &graph) == 3);
  TASSERT(kmer_get_covg("TAGAACGTTCCCTACACGT", &graph) == 3);

  // Update statistics
  graph_info_update_stats(&graph.ginfo[0], &stats);

  double mean_read_length = ((double)total_seq/contigs_loaded)+0.5;

  size_t g_total_seq = graph.ginfo[0].total_sequence;
  size_t g_mean_read_length = graph.ginfo[0].mean_read_length;

  // test_status("g_total_seq: %zu total_seq: %zu", g_total_seq, total_seq);
  // test_status("g_mean_read_length: %zu mean_read_length: %f",
  //             g_mean_read_length, mean_read_length);

  TASSERT2(g_total_seq == total_seq, "%zu %zu", g_total_seq, total_seq);
  TASSERT(fabs(g_mean_read_length - mean_read_length) <= 0.5);

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  db_graph_dealloc(&graph);
}
