#include "global.h"
#include "all_tests.h"
#include "build_graph.h"
#include "correct_alignment.h"
#include "generate_paths.h"
#include "db_alignment.h"

static void _check_correct_aln(char *seq1, char *seq2,
                               char **ans, size_t num_ans,
                               dBAlignment *aln,
                               CorrectAlnWorker *corrector,
                               const CorrectAlnParam params,
                               const dBGraph *graph, StrBuf *sbuf)
{
  size_t i, seqlen = strlen(seq1);
  dBNodeBuffer *nbuf;

  // Fake read
  char empty[10] = "", rname[20] = "Example";
  read_t r1 = {.name = {.b = rname, .end = strlen(rname), .size = 10},
               .seq  = {.b = seq1, .end = seqlen, .size = seqlen+1},
               .qual = {.b = empty, .end = 0, .size = 1}};

  read_t r2 = {.name = {.b = rname, .end = strlen(rname), .size = 10},
               .seq  = {.b = empty, .end = 0, .size = 1},
               .qual = {.b = empty, .end = 0, .size = 1}};
  read_t *r2ptr = NULL;

  if(seq2) {
    r2.seq.b = seq2;
    r2.seq.end = strlen(seq2);
    r2.seq.size = r2.seq.end+1;
    r2ptr = &r2;
  }

  db_alignment_from_reads(aln, &r1, r2ptr, 0, 0, 0, graph, -1);
  correct_alignment_init(corrector, aln, params);

  for(i = 0; i < num_ans; i++)
  {
    nbuf = correct_alignment_nxt(corrector);
    TASSERT(nbuf != NULL);

    strbuf_ensure_capacity(sbuf, nbuf->len+MAX_KMER_SIZE+1);
    db_nodes_to_str(nbuf->data, nbuf->len, graph, sbuf->buff);

    TASSERT2(strcmp(sbuf->buff, ans[i]) == 0,
             "%zu) exp: %s got: %s", i, ans[i], sbuf->buff);
  }

  // Next alignment should be NULL
  nbuf = correct_alignment_nxt(corrector);
  TASSERT(nbuf == NULL);
}


static void test_correct_aln_no_paths()
{
  test_status("[CorrectAln] Testing correct_aln with no paths...");

  // mutations:                            **                 *
  char seq[] = "ATGCATGTTGACCAAATAAGTCAC""TGTGGGAGCCACGTAAAGCGTTCGCACCGATTTGTG";
  char mu0[] =     "ATGTTGACCAAATAAGTCAC""TGTCCGAGCCACGTAAAGCGTTCGCACC";
  char re0[] =     "ATGTTGACCAAATAAGTCAC""TGTGGGAGCCACGTAAAGCGTTCGCACC";
  //                                    v                       *
  char mu1[] =     "ATGTTGACCAAATAAGTCA" "TGTGGGAGCCACGTAAAGCGTTAGCACCGATTTGTG";
  char re1[] =     "ATGTTGACCAAATAAGTCAC""TGTGGGAGCCACGTAAAGCGTTCGCACCGATTTGTG";

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  // Start alignment
  dBAlignment aln;
  CorrectAlnWorker corrector;
  StrBuf sbuf;

  // Set up alignment correction params
  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .ins_gap_min = 0, .ins_gap_max = 0,
                            .one_way_gap_traverse = true, .max_context = 10,
                            .gap_variance = 0.1, .gap_wiggle = 5};

  char *gseqs[1] = {seq}, *alns[1];

  // Construct graph and paths
  _construct_graph_with_paths(&graph, kmer_size, ncols, gseqs, 1, params);

  TASSERT(graph.pstore.num_of_paths == 0);

  db_alignment_alloc(&aln);
  correct_aln_worker_alloc(&corrector, &graph);
  strbuf_alloc(&sbuf, 1024);

  alns[0] = re0;
  _check_correct_aln(mu0, NULL, alns, 1, &aln, &corrector, params, &graph, &sbuf);

  alns[0] = re1;
  _check_correct_aln(mu1, NULL, alns, 1, &aln, &corrector, params, &graph, &sbuf);

  strbuf_dealloc(&sbuf);
  correct_aln_worker_dealloc(&corrector);
  db_alignment_dealloc(&aln);
  _deconstruct_graph_with_paths(&graph);
}

static void test_contig_ends_agree()
{
  test_status("[CorrectAln] Testing correct_aln with contig end check...");

  //    1         2         3 < mutations
  // a --+      +-a-+      +YY a
  //      \    /     \    Y    
  //       +X-+       +-Y+     
  //      X    \     /    \    
  // b XX+      +-b-+      +-- b

  // Read pair XXX <gap> YYYY
  // Should fail graph_walker_agrees_contig() test

  // mutations:   1               2                       3
  char seqa[] = "CCGATTAAAGGGTTACTATAGCACAGGAATGGTCTGGCCTGTAAGAAGTCCAGCTTC"; // a
  char seqb[] = "CAGATTAAAGGGTTACTGTAGCACAGGAATGGTCTGGCCTGTAAGATGTCCAGCTTC"; // b

  // read crosses from b -> a
  char r1a[] =  "CCGATTAAAGGGTT";
  char r1b[] =  "CAGATTAAAGGGTT";
  char r2a[] =                                     "GGCCTGTAAGAAGTCCAGCTTC";
  char r2b[] =                                     "GGCCTGTAAGATGTCCAGCTTC";

// ATTGAAGCTGG

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  // Start alignment
  dBAlignment aln;
  CorrectAlnWorker corrector;
  StrBuf sbuf;

  // Set up alignment correction params
  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .ins_gap_min = 0, .ins_gap_max = 30,
                            .one_way_gap_traverse = true, .max_context = 10,
                            .gap_variance = 0.1, .gap_wiggle = 5};

  char *seqs[2] = {seqa, seqb}, *alns[2];

  // Construct graph and paths
  _construct_graph_with_paths(&graph, kmer_size, ncols, seqs, 2, params);

  db_alignment_alloc(&aln);
  correct_aln_worker_alloc(&corrector, &graph);
  strbuf_alloc(&sbuf, 1024);

  // Check number of kmers in the graph
  size_t expnkmers = strlen(seqa)+1-kmer_size + 2+11+11;
  TASSERT2(graph.ht.num_kmers == expnkmers,
           "%zu vs %zu", (size_t)graph.ht.num_kmers, expnkmers);

  // Check number of paths
  TASSERT(graph.pstore.num_of_paths == 8);
  TASSERT(graph.pstore.num_kmers_with_paths == 8);
  TASSERT(graph.pstore.num_col_paths == 8);

  size_t t;

  for(t = 0; t < 2; t++)
  {
    params.one_way_gap_traverse = (t == 0);

    correct_aln_stats_zero(&corrector.gapstats);
    alns[0] = seqa;
    _check_correct_aln(r1a, r2a, alns, 1, &aln, &corrector, params, &graph, &sbuf);
    TASSERT(corrector.gapstats.num_gap_successes == 1);

    correct_aln_stats_zero(&corrector.gapstats);
    alns[0] = seqb;
    _check_correct_aln(r1b, r2b, alns, 1, &aln, &corrector, params, &graph, &sbuf);
    TASSERT(corrector.gapstats.num_gap_successes == 1);

    correct_aln_stats_zero(&corrector.gapstats);
    alns[0] = r1a;
    alns[1] = r2b;
    _check_correct_aln(r1a, r2b, alns, 2, &aln, &corrector, params, &graph, &sbuf);
    TASSERT(corrector.gapstats.num_gap_successes == 0);
    TASSERT(corrector.gapstats.num_gaps_disagreed > 0);

    correct_aln_stats_zero(&corrector.gapstats);
    alns[0] = r1b;
    alns[1] = r2a;
    _check_correct_aln(r1b, r2a, alns, 2, &aln, &corrector, params, &graph, &sbuf);
    TASSERT(corrector.gapstats.num_gap_successes == 0);
    TASSERT(corrector.gapstats.num_gaps_disagreed > 0);
  }

  strbuf_dealloc(&sbuf);
  correct_aln_worker_dealloc(&corrector);
  db_alignment_dealloc(&aln);
  _deconstruct_graph_with_paths(&graph);
}

void test_corrected_aln()
{
  test_status("[CorrectAln] Testing correct_alignment.c");
  test_correct_aln_no_paths();
  test_contig_ends_agree();
}
