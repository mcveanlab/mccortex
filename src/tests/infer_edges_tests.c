#include "global.h"
#include "all_tests.h"

#include "db_node.h"
#include "db_graph.h"
#include "infer_edges.h"

#define node_in_col(g,n,c) \
        ((n).key != HASH_NOT_FOUND && db_node_has_col((g), (n).key, (c)))

#define LEFTNODE 8
#define RGHTNODE 4
#define LEFTEDGE 2
#define RGHTEDGE 1

#define ALLEDGES (LEFTNODE|RGHTNODE|LEFTEDGE|RGHTEDGE)

// Requires string of length kmer_size + 1
// Returns:
//   0b8421
// 8 - has left node
// 4 - has right node
// 2 - has edge from left node
// 1 - has edge from right node
static uint8_t get_edges(const char *seq, const dBGraph *db_graph, size_t col)
{
  ctx_assert(strlen(seq) == db_graph->kmer_size+1);
  dBNode n0 = db_graph_find_str(db_graph, seq);
  dBNode n1 = db_graph_find_str(db_graph, seq+1);
  uint8_t key = (node_in_col(db_graph,n0,col) ? 8 : 0) |
                (node_in_col(db_graph,n1,col) ? 4 : 0);
  if(key != 12) return key; // return if we don't have both nodes
  Edges edges0 = db_node_edges(db_graph, n0.key, col);
  Edges edges1 = db_node_edges(db_graph, n1.key, col);
  Nucleotide lhs_nuc = dna_char_to_nuc(seq[0]);
  Nucleotide rhs_nuc = dna_char_to_nuc(seq[db_graph->kmer_size]);
  uint8_t e0 = edges_has_edge(edges0, rhs_nuc, n0.orient) ? 2 : 0;
  uint8_t e1 = edges_has_edge(edges1, dna_nuc_complement(lhs_nuc), !n1.orient) ? 1 : 0;
  return (key | e0 | e1);
}

static void simple_test()
{
  // Construct 4 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 5;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 2000);
  // Graph data
  graph.bktlocks = ctx_calloc(roundup_bits2bytes(graph.ht.num_of_buckets),
                              sizeof(graph.bktlocks[0]));
  graph.col_edges = ctx_calloc(graph.ht.capacity * ncols, sizeof(Edges));
  graph.node_in_cols = ctx_calloc(roundup_bits2bytes(graph.ht.capacity)*ncols,
                                  sizeof(graph.node_in_cols[0]));

  // TAACAATGACT -> AACAATGACTC -> ACAATGACTCC
  //                            -> ACAATGACTCG
  //
  //  - => edge, * => missing edge
  //
  // Sample 0 has all nodes, all edges
  // 0-0-0
  //    -0
  // Sample 1 has all nodes, no edges
  // 1*1*1  ->  1-1-1
  //    *1         -1
  // Sample 2 has some nodes, missing one edge
  // 2-2 .  ->  2-2 .
  //    *2         -2
  // Sample 3 has some nodes, missing no edges
  // 3 . 3  ->  3 . 3
  //     3          3
  // Samle 4 has all nodes, missing two edges
  // 4*4*4  -> 4-4-4
  //    -4        -4

  //                   V           *
  const char seq0[] = "TAACAATGACTCC";
  const char seq1[] =  "AACAATGACTCG";

  // Add random kmer that no one has edges to/from
  const char rnd[] = "GGACTTTTTAA";
  build_graph_from_str_mt(&graph, 0, rnd, 11);
  build_graph_from_str_mt(&graph, 1, rnd, 11);
  build_graph_from_str_mt(&graph, 2, rnd, 11);
  build_graph_from_str_mt(&graph, 3, rnd, 11);

  // Sample 0
  build_graph_from_str_mt(&graph, 0, seq0, 13);
  build_graph_from_str_mt(&graph, 0, seq1, 12);

  // Sample 1
  build_graph_from_str_mt(&graph, 1, seq0,   11);
  build_graph_from_str_mt(&graph, 1, seq0+1, 11);
  build_graph_from_str_mt(&graph, 1, seq0+2, 11);
  build_graph_from_str_mt(&graph, 1, seq1+1, 11);

  // Sample 2
  build_graph_from_str_mt(&graph, 2, seq0,   12);
  build_graph_from_str_mt(&graph, 2, seq1+1, 11);

  // Sample 3
  build_graph_from_str_mt(&graph, 3, seq0,   11);
  build_graph_from_str_mt(&graph, 3, seq0+2, 11);
  build_graph_from_str_mt(&graph, 3, seq1+1, 11);

  // Sample 4
  build_graph_from_str_mt(&graph, 4, seq0,   11);
  build_graph_from_str_mt(&graph, 4, seq0+1, 11);
  build_graph_from_str_mt(&graph, 4, seq0+2, 11);
  build_graph_from_str_mt(&graph, 4, seq1,   12);

  TASSERT(graph.ht.num_kmers == 5);

  // First edge
  TASSERT(get_edges("TAACAATGACTC", &graph, 0) == ALLEDGES);
  TASSERT(get_edges("TAACAATGACTC", &graph, 1) == (LEFTNODE|RGHTNODE));
  TASSERT(get_edges("TAACAATGACTC", &graph, 2) == ALLEDGES);
  TASSERT(get_edges("TAACAATGACTC", &graph, 3) == LEFTNODE);
  TASSERT(get_edges("TAACAATGACTC", &graph, 4) == (LEFTNODE|RGHTNODE));

  // Second edge
  TASSERT(get_edges("AACAATGACTCC", &graph, 0) == ALLEDGES);
  TASSERT(get_edges("AACAATGACTCC", &graph, 1) == (LEFTNODE|RGHTNODE));
  TASSERT(get_edges("AACAATGACTCC", &graph, 2) == LEFTNODE);
  TASSERT(get_edges("AACAATGACTCC", &graph, 3) == RGHTNODE);
  TASSERT(get_edges("AACAATGACTCC", &graph, 4) == (LEFTNODE|RGHTNODE));

  // last edge
  TASSERT(get_edges("AACAATGACTCG", &graph, 0) == ALLEDGES);
  TASSERT(get_edges("AACAATGACTCG", &graph, 1) == (LEFTNODE|RGHTNODE));
  TASSERT(get_edges("AACAATGACTCG", &graph, 2) == (LEFTNODE|RGHTNODE));
  TASSERT(get_edges("AACAATGACTCG", &graph, 3) == RGHTNODE);
  TASSERT(get_edges("AACAATGACTCG", &graph, 4) == ALLEDGES);

  size_t num_modified_nodes;
  num_modified_nodes = infer_edges(2, true, &graph);

  // All four nodes modified
  TASSERT2(num_modified_nodes == 4, "modified: %zu", num_modified_nodes);

  // Check edges have been added: (LEFTNODE|RGHTNODE) => ALLEDGES
  // First edge
  TASSERT(get_edges("TAACAATGACTC", &graph, 0) == ALLEDGES);
  TASSERT(get_edges("TAACAATGACTC", &graph, 1) == ALLEDGES);
  TASSERT(get_edges("TAACAATGACTC", &graph, 2) == ALLEDGES);
  TASSERT(get_edges("TAACAATGACTC", &graph, 3) == LEFTNODE);
  TASSERT(get_edges("TAACAATGACTC", &graph, 4) == ALLEDGES);

  // Second edge
  TASSERT(get_edges("AACAATGACTCC", &graph, 0) == ALLEDGES);
  TASSERT(get_edges("AACAATGACTCC", &graph, 1) == ALLEDGES);
  TASSERT(get_edges("AACAATGACTCC", &graph, 2) == LEFTNODE);
  TASSERT(get_edges("AACAATGACTCC", &graph, 3) == RGHTNODE);
  TASSERT(get_edges("AACAATGACTCC", &graph, 4) == ALLEDGES);

  // last edge
  TASSERT(get_edges("AACAATGACTCG", &graph, 0) == ALLEDGES);
  TASSERT(get_edges("AACAATGACTCG", &graph, 1) == ALLEDGES);
  TASSERT(get_edges("AACAATGACTCG", &graph, 2) == ALLEDGES);
  TASSERT(get_edges("AACAATGACTCG", &graph, 3) == RGHTNODE);
  TASSERT(get_edges("AACAATGACTCG", &graph, 4) == ALLEDGES);

  ctx_free(graph.bktlocks);
  ctx_free(graph.col_edges);
  ctx_free(graph.node_in_cols);
  db_graph_dealloc(&graph);
}

void test_infer_edges_tests()
{
  test_status("Testing infer_edges...");
  simple_test();
}
