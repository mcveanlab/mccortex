#include "global.h"
#include "all_tests.h"
#include "binary_kmer.h"
#include "db_node.h"
#include "db_unitig.h"
#include "build_graph.h"

#include "bit_array/bit_macros.h"

#define UNODEBUF 200

static void unitig_from_kmer(hkey_t hkey, dBNodeBuffer *nbuf,
                             uint64_t *visited, const dBGraph *graph,
                             const char **ans, size_t n)
{
  size_t i;
  char tmpstr[UNODEBUF];

  if(!bitset_get(visited, hkey))
  {
    db_node_buf_reset(nbuf);
    db_unitig_fetch(hkey, nbuf, graph);
    for(i = 0; i < nbuf->len; i++) bitset_set(visited, nbuf->b[i].key);

    db_unitig_normalise(nbuf->b, nbuf->len, graph);

    TASSERT(nbuf->len < UNODEBUF);
    db_nodes_to_str(nbuf->b, nbuf->len, graph, tmpstr);
    for(i = 0; i < n && strcmp(tmpstr,ans[i]) != 0; i++);

    TASSERT2(i < n, "Got: %s", tmpstr);
  }
}

static void pull_out_unitigs(const char **seq, const char **ans, size_t n,
                             const dBGraph *graph)
{
  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 1024);

  // 1. Check pulling out unitigs works for iterating over the graph
  uint64_t *visited;
  visited = ctx_calloc(roundup_bits2words64(graph->ht.capacity), 8);
  HASH_ITERATE(&graph->ht, unitig_from_kmer,
               &nbuf, visited, graph, ans, n);
  ctx_free(visited);

  // 2. Check pulling out unitigs works when we iterate over inputs
  size_t i, j, len;
  dBNode node;
  char tmpstr[UNODEBUF];

  for(i = 0; i < n; i++) {
    len = strlen(seq[i]);
    for(j = 0; j+graph->kmer_size <= len; j++)
    {
      // Find node
      node = db_graph_find_str(graph, seq[i]+j);
      TASSERT(node.key != HASH_NOT_FOUND);

      // Fetch unitig
      db_node_buf_reset(&nbuf);
      db_unitig_fetch(node.key, &nbuf, graph);
      db_unitig_normalise(nbuf.b, nbuf.len, graph);

      // Compare
      TASSERT(nbuf.len < UNODEBUF);
      db_nodes_to_str(nbuf.b, nbuf.len, graph, tmpstr);
      if(strcmp(tmpstr, ans[i]) != 0) {
        test_status("Got: %s from ans[i]:%s\n", tmpstr, ans[i]);
      }
      TASSERT(strcmp(tmpstr, ans[i]) == 0);
    }
  }

  db_node_buf_dealloc(&nbuf);
}

void test_db_unitig()
{
  test_status("testing db_unitig_fetch()...");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 19, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 1024,
                 DBG_ALLOC_EDGES | DBG_ALLOC_COVGS | DBG_ALLOC_BKTLOCKS);

  #define NSEQ 7

  const char *seq[NSEQ]
   = {"AGAGAGAGAGAGAGAGAGAGAGAG",
      "AAAAAAAAAAAAAAAAAAAAAAAAAA",
      "ATATATATATATATATATATATATATAT",
      "CGTTCGCGCATGGCCCACG",
      "GAACCAATCGGTCGACTGT",
      "CCCCGCAAAGTCCACTTAGTGTAAGGTACAAATTCTGCAGAGTTGCTGGATCAGCGATAC",
      "TCAATCCGATAGCAACCCGGTCCAA""TCAATCCGATAGCAACCCGGTCCAA"};

  const char *ans[NSEQ]
   = {"AGAGAGAGAGAGAGAGAGAG", // key AGAGAGAGAGAGAGAGAGA < CTCTCTCTCTCTCTCTCTC
      "AAAAAAAAAAAAAAAAAAA",
      "ATATATATATATATATATA",
      "CGTGGGCCATGCGCGAACG",
      "ACAGTCGACCGATTGGTTC",
      "CCCCGCAAAGTCCACTTAGTGTAAGGTACAAATTCTGCAGAGTTGCTGGATCAGCGATAC",
      "AACCCGGTCCAATCAATCCGATAGCAACCCGGTCCAATCAATC"};

  // Load all seq into colour 0
  size_t i;
  for(i = 0; i < NSEQ; i++)
    build_graph_from_str_mt(&graph, 0, seq[i], strlen(seq[i]), false);

  pull_out_unitigs(seq, ans, NSEQ, &graph);

  db_graph_dealloc(&graph);
}
