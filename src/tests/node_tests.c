#include "global.h"
#include "all_tests.h"
#include <ctype.h>

#include "util.h"
#include "db_graph.h"
#include "db_node.h"
#include "build_graph.h"

static void edge_check(hkey_t hkey, const dBGraph *db_graph, size_t col)
{
  const BinaryKmer bkmer = db_node_get_bkmer(db_graph, hkey);
  const Edges edges = db_node_get_edges(db_graph, hkey, col);

  dBNode nodes[4];
  Nucleotide nucs[4];
  size_t i, n, or;

  for(or = 0; or < 2; or++) {
    Edges e = 0;
    n = db_graph_next_nodes(db_graph, bkmer, or, edges, nodes, nucs);
    for(i = 0; i < n; i++) e |= nuc_orient_to_edge(nucs[i], or);
    TASSERT(edges_with_orientation(e,or) == edges_with_orientation(edges,or));
  }
}

static void test_db_graph_next_nodes()
{
  test_status("Testing db_graph_next_nodes() vs edges_ functions");

  // Construct 2 colour graph with kmer-size=11
  dBGraph graph;
  size_t col, kmer_size = 11, ncols = 2;
  char seq[60];

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 1024,
                 DBG_ALLOC_EDGES | DBG_ALLOC_COVGS | DBG_ALLOC_BKTLOCKS);

  // Copy a random and shared piece of sequence to both colours
  for(col = 0; col < 2; col++) {
    dna_rand_str(seq, 59);
    build_graph_from_str_mt(&graph, col, seq, strlen(seq), false);
    strcpy(seq, "CTTTCTTATCTGGAACCAGCTTTGCGGGGATGGAGTGTAACCTTGACAATGGGTCCTGC");
    build_graph_from_str_mt(&graph, col, seq, strlen(seq), false);
  }

  HASH_ITERATE(&graph.ht, edge_check, &graph, 0);
  HASH_ITERATE(&graph.ht, edge_check, &graph, 1);

  db_graph_dealloc(&graph);
}

#define MAXLEN 300
#define NLOOP 300

static void test_left_shift()
{
  test_status("Testing db_nodes_left_shift()");

  dBNode nodes0[MAXLEN], nodes1[MAXLEN];

  size_t i, j, len, shift;
  for(i = 0; i < NLOOP; i++)
  {
    len = rand() % MAXLEN;
    shift = len ? rand() % len : 0;

    // Generate rand array
    for(j = 0; j < len; j++) {
      nodes0[j].key = rand();
      nodes0[j].orient = (nodes0[j].key>>3)&1;
    }
    // Make two copies
    memcpy(nodes1, nodes0, len * sizeof(dBNode));

    // Use slow shift
    db_nodes_reverse(nodes0, shift);
    db_nodes_reverse(nodes0+shift, len-shift);
    db_nodes_reverse(nodes0, len);

    // Use fast shift
    db_nodes_left_shift(nodes1, len, shift);

    // Compare
    for(j = 0; j < len; j++)
      TASSERT(db_nodes_are_equal(nodes0[j], nodes1[j]));
  }
}

void test_db_node()
{
  test_db_graph_next_nodes();
  test_left_shift();
}
