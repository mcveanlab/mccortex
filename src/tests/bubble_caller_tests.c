#include "global.h"
#include "all_tests.h"

#include "bubble_caller.h"
#include "genotyping.h" // Tested here for now

static void _check_alleles(GraphCache *cache, GCacheStepPtrBuf *steps,
                           const char **alleles, size_t num_alleles,
                           dBNodeBuffer *nbuf, StrBuf *sbuf)
{
  TASSERT2(steps->len == num_alleles, "Number of alleles doesn't match");

  size_t i, j;
  for(i = 0; i < steps->len; i++)
  {
    db_node_buf_reset(nbuf);
    gc_step_fetch_nodes(cache, steps->b[i], nbuf);
    strbuf_ensure_capacity(sbuf, nbuf->len+MAX_KMER_SIZE+1);
    db_nodes_to_str(nbuf->b, nbuf->len, cache->db_graph, sbuf->b);

    // Find this node
    for(j = 0; j < num_alleles && strcasecmp(sbuf->b,alleles[j]); j++) {}
    TASSERT2(j < num_alleles, "Couldn't find allele: %s", sbuf->b);
  }
}

// `nbuf` and `sbuf` are temporary variables used by this function
static void _call_bubble(BubbleCaller *caller,
                         const char *flank5p, const char *flank3p,
                         const char **alleles, size_t num_alleles,
                         dBNodeBuffer *nbuf, StrBuf *sbuf)
{
  const dBGraph *graph = caller->db_graph;
  const size_t kmer_size = graph->kmer_size;

  dBNode node5p = db_graph_find_str(graph, flank5p+strlen(flank5p)-kmer_size);
  dBNode node3p = db_graph_find_str(graph, flank3p);
  TASSERT(node5p.key != HASH_NOT_FOUND);
  TASSERT(node3p.key != HASH_NOT_FOUND);

  Edges edges5p = db_node_get_edges_union(graph, node5p.key);
  Edges edges3p = db_node_get_edges_union(graph, node3p.key);
  TASSERT(edges_get_outdegree(edges5p, node5p.orient) > 1);
  TASSERT(edges_get_indegree(edges3p, node3p.orient) > 1);

  find_bubbles(caller, node5p);

  GCacheUnitig *snode3p;
  Orientation snorient3p;
  GCacheStepPtrBuf *stepbuf;

  // Get 3p flank and orientation
  snode3p = graph_cache_find_unitig(&caller->cache, node3p);
  TASSERT(snode3p != NULL);
  snorient3p = gc_unitig_get_orient(&caller->cache, snode3p, node3p);

  find_bubbles_ending_with(caller, snode3p);

  stepbuf = (snorient3p == FORWARD ? &caller->spp_forward : &caller->spp_reverse);

  _check_alleles(&caller->cache, stepbuf, alleles, num_alleles, nbuf, sbuf);
}

// Load each sequence into a separate colour
static void test_bubbles(dBGraph *graph, const char **seqs, size_t nseqs,
                         const char *flank5p, const char *flank3p,
                         const char **alleles, size_t nalleles)
{
  db_graph_reset(graph);

  TASSERT(graph->num_of_cols >= nseqs);

  size_t i;
  for(i = 0; i < nseqs; i++)
    build_graph_from_str_mt(graph, i, seqs[i], strlen(seqs[i]));

  graph->num_of_cols_used = MAX2(graph->num_of_cols_used, 1);

  StrBuf sbuf;
  dBNodeBuffer nbuf;
  strbuf_alloc(&sbuf, 128);
  db_node_buf_alloc(&nbuf, 128);

  BubbleCallingPrefs prefs = {.max_allele_len = 100, .max_flank_len = 100,
                              .haploid_cols = NULL, .nhaploid_cols = 0,
                              .remove_serial_bubbles = true};

  BubbleCaller *caller = bubble_callers_new(1, &prefs, NULL, graph);

  _call_bubble(caller, flank5p, flank3p, alleles, nalleles, &nbuf, &sbuf);

  strbuf_dealloc(&sbuf);
  db_node_buf_dealloc(&nbuf);
  bubble_callers_destroy(caller, 1);
}

void test_bubble_caller()
{
  test_status("Testing bubble calling...");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  const size_t kmer_size = 11, ncols = 3;

  // Create graph
  db_graph_alloc(&graph, kmer_size, ncols, 1, 2000,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL | DBG_ALLOC_BKTLOCKS);

  //   mutations:                      x
  const char *seqs0[] = {"AGGGATAAAACTCTGTACTGGATCTCCCT",
                         "AGGGATAAAACTCTcTACTGGATCTCCCT"};
  const char flank5p0[] = "AGGGATAAAACTCT";
  const char flank3p0[] = "TACTGGATCTCCCT";
  const char *alleles0[] = {"ATAAAACTCTGTACTGGATCT", "ATAAAACTCTcTACTGGATCT"};

  test_bubbles(&graph, seqs0, 2, flank5p0, flank3p0, alleles0, 2);

  //   mutations:                     x                  y
  const char *seqs1[] = {"CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA",
                         "CCCGTAGGTAAGtGCGTTAGTGCAAGGCCACATTGGGACACGAGTTGATA",
                         "CCCGTAGGTAAGGGCGTTAGTGCAAGGCCACtTTGGGACACGAGTTGATA"};

  // forwards
  const char flank5p1a[] = "CCCGTAGGTAAG";
  const char flank3p1a[] = "GCGTTAGTGCAAGGCCAC";
  const char *alleles1a[] = {"CGTAGGTAAGGGCGTTAGTGC", "CGTAGGTAAGtGCGTTAGTGC"};

  const char flank5p1b[] = "GCGTTAGTGCAAGGCCAC";
  const char flank3p1b[] = "TTGGGACACGAGTTGATA";
  const char *alleles1b[] = {"GCAAGGCCACATTGGGACACG", "GCAAGGCCACtTTGGGACACG"};

  test_bubbles(&graph, seqs1, 3, flank5p1a, flank3p1a, alleles1a, 2);
  test_bubbles(&graph, seqs1, 3, flank5p1b, flank3p1b, alleles1b, 2);

  // reverse
  // mutations:        y                  x
  // TATCAACTCGTGTCCCAATGTGGCCTTGCACTAACGCCCTTACCTACGGG
  // TATCAACTCGTGTCCCAATGTGGCCTTGCACTAACGCaCTTACCTACGGG
  // TATCAACTCGTGTCCCAAaGTGGCCTTGCACTAACGCCCTTACCTACGGG
  //
  const char flank5p1c[] = "GTGGCCTTGCACTAACGC";
  const char flank3p1c[] = "CTTACCTACGGG";
  const char *alleles1c[] = {"GCACTAACGCCCTTACCTACG", "GCACTAACGCaCTTACCTACG"};

  const char flank5p1d[] = "TATCAACTCGTGTCCCAA";
  const char flank3p1d[] = "GTGGCCTTGCACTAACGC";
  const char *alleles1d[] = {"CGTGTCCCAATGTGGCCTTGC", "CGTGTCCCAAaGTGGCCTTGC"};

  test_bubbles(&graph, seqs1, 3, flank5p1c, flank3p1c, alleles1c, 2);
  test_bubbles(&graph, seqs1, 3, flank5p1d, flank3p1d, alleles1d, 2);

  db_graph_dealloc(&graph);
}
