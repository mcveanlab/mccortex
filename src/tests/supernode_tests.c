#include "global.h"
#include "all_tests.h"
#include "binary_kmer.h"
#include "db_node.h"
#include "supernode.h"
#include "build_graph.h"

#include "bit_array/bit_macros.h"

#define SNODEBUF 200

static void supernode_from_kmer(hkey_t hkey, dBNodeBuffer *nbuf,
                                uint64_t *visited, const dBGraph *graph,
                                const char **ans, size_t n)
{
  size_t i;
  char tmpstr[SNODEBUF];

  if(!bitset_get(visited, hkey))
  {
    db_node_buf_reset(nbuf);
    supernode_find(hkey, nbuf, graph);
    for(i = 0; i < nbuf->len; i++) bitset_set(visited, nbuf->data[i].key);

    supernode_normalise(nbuf->data, nbuf->len, graph);

    TASSERT(nbuf->len < SNODEBUF);
    db_nodes_to_str(nbuf->data, nbuf->len, graph, tmpstr);
    for(i = 0; i < n && strcmp(tmpstr,ans[i]) != 0; i++);

    TASSERT2(i < n, "Got: %s", tmpstr);
  }
}

static void pull_out_supernodes(const char **seq, const char **ans, size_t n,
                                const dBGraph *graph)
{
  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 1024);

  // 1. Check pulling out supernodes works for iterating over the graph
  uint64_t *visited;
  visited = ctx_calloc(roundup_bits2words64(graph->ht.capacity), 8);
  HASH_ITERATE(&graph->ht, supernode_from_kmer,
               &nbuf, visited, graph, ans, n);
  ctx_free(visited);

  // 2. Check pulling out supernodes works when we iterate over inputs
  size_t i, j, len;
  BinaryKmer bkmer;
  dBNode node;
  char tmpstr[SNODEBUF];

  for(i = 0; i < n; i++) {
    len = strlen(seq[i]);
    for(j = 0; j+graph->kmer_size <= len; j++)
    {
      // Find node
      bkmer = binary_kmer_from_str(seq[i]+j, graph->kmer_size);
      node = db_graph_find(graph, bkmer);
      TASSERT(node.key != HASH_NOT_FOUND);

      // Fetch supernode
      db_node_buf_reset(&nbuf);
      supernode_find(node.key, &nbuf, graph);
      supernode_normalise(nbuf.data, nbuf.len, graph);

      // Compare
      TASSERT(nbuf.len < SNODEBUF);
      db_nodes_to_str(nbuf.data, nbuf.len, graph, tmpstr);
      if(strcmp(tmpstr, ans[i]) != 0) {
        test_status("Got: %s from ans[i]:%s\n", tmpstr, ans[i]);
      }
      TASSERT(strcmp(tmpstr, ans[i]) == 0);
    }
  }

  db_node_buf_dealloc(&nbuf);
}

void test_supernode()
{
  test_status("testing supernode_find()...");

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
    build_graph_from_str_mt(&graph, 0, seq[i], strlen(seq[i]));

  pull_out_supernodes(seq, ans, NSEQ, &graph);

  db_graph_dealloc(&graph);
}
