#include "global.h"
#include "all_tests.h"

#include "db_node.h"
#include "db_graph.h"
#include "graph_crawler.h"


void test_graph_crawler()
{
  test_status("Testing graph crawler...");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  const size_t kmer_size = 11, ncols = 3;

  db_graph_alloc(&graph, kmer_size, ncols, 1, 2048);
  // Graph data
  graph.bktlocks = ctx_calloc(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = ctx_calloc(graph.ht.capacity, sizeof(Edges));
  graph.node_in_cols = ctx_calloc(roundup_bits2words64(graph.ht.capacity)*ncols, 8);

  char graphseq[3][77] =
//           <               X                 X              X...............
{"GTTCCAGAGCGGAGGTCTCCCAACAACATGGTATAAGTTGTCTAGCCCCGGTTCGCGCGGGTACTTCTTACAGCGC",
 "GTTCCAGAGCGGAGGTCTCCCAACAACTTGGTATAAGTTGTCTAGTCCCGGTTCGCGCGGCATTTCAGCATTGTTA",
 "GTTCCAGAGCGCGACAGAGTGCATATCACGCTAAGCACAGCCCTCTTCTATCTGCTTTTAAATGGATCAATAATCG"};

  build_graph_from_str_mt(&graph, 0, graphseq[0], strlen(graphseq[0]));
  build_graph_from_str_mt(&graph, 1, graphseq[1], strlen(graphseq[1]));
  build_graph_from_str_mt(&graph, 2, graphseq[2], strlen(graphseq[2]));

  // Crawl graph
  GraphCrawler crawler;
  graph_crawler_alloc(&crawler, &graph);

  dBNode node = db_graph_find_str(&graph, graphseq[0]);
  dBNode next_node = db_graph_find_str(&graph, graphseq[0]+1);
  TASSERT(node.key != HASH_NOT_FOUND);
  TASSERT(next_node.key != HASH_NOT_FOUND);

  BinaryKmer bkey = db_node_get_bkmer(&graph, node.key);
  Edges edges = db_node_get_edges(&graph, node.key, 0);

  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  size_t i, p, num_next, next_idx;

  num_next = db_graph_next_nodes(&graph, bkey, node.orient, edges,
                                 next_nodes, next_nucs);

  next_idx = 0;
  while(next_idx < num_next && !db_nodes_are_equal(next_nodes[next_idx],next_node))
    next_idx++;

  TASSERT(next_idx < num_next && db_nodes_are_equal(next_nodes[next_idx],next_node));

  // Crawl in all colours
  graph_crawler_fetch(&crawler, node, next_nodes, next_nucs, next_idx, num_next,
                      NULL, graph.num_of_cols, NULL, NULL, NULL);

  TASSERT2(crawler.num_paths == 2, "crawler.num_paths: %u", crawler.num_paths);

  // Fetch paths
  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 16);
  StrBuf sbuf;
  strbuf_alloc(&sbuf, 128);

  for(p = 0; p < crawler.num_paths; p++) {
    db_node_buf_reset(&nbuf);
    graph_crawler_get_path_nodes(&crawler, p, &nbuf);
    strbuf_ensure_capacity(&sbuf, nbuf.len+graph.kmer_size);
    sbuf.end = db_nodes_to_str(nbuf.data, nbuf.len, &graph, sbuf.b);
    for(i = 0; i < 3 && strcmp(graphseq[i]+1,sbuf.b) != 0; i++) {}
    TASSERT2(i < 3, "seq: %s", sbuf.b);
    TASSERT2(sbuf.end == 75, "sbuf.end: %zu", sbuf.end);
    TASSERT2(nbuf.len == 65, "nbuf.len: %zu", nbuf.len);
  }

  strbuf_dealloc(&sbuf);
  db_node_buf_dealloc(&nbuf);

  graph_crawler_dealloc(&crawler);

  db_graph_dealloc(&graph);
}
