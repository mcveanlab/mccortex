#include "global.h"
#include "all_tests.h"
#include "sorted_path_set.h"
#include "binary_seq.h"

#include "string_buffer.h"

//
// Faked sorted path set
//

static void _fake_sorted_path_set(SortedPathSet *set, size_t cbytes,
                                  const char **paths, size_t npaths)
{
  qsort(paths, npaths, sizeof(char*), cmp_charptr);

  set->hkey = 0; // this is not important
  set->cbytes = cbytes;
  bytebuf_reset(&set->seqs);
  pentrybuf_reset(&set->members);

  size_t i, plen, total_bytes = 0;
  uint8_t *seq;

  for(i = 0; i < npaths; i++) total_bytes += cbytes + strlen(paths[i]);

  bytebuf_ensure_capacity(&set->seqs, total_bytes);
  pentrybuf_ensure_capacity(&set->members, npaths);

  for(i = 0; i < npaths; i++)
  {
    plen = strlen(paths[i]);
    seq = set->seqs.data + set->seqs.len;
    seq[0] = 1; // single colour set in colset
    PathEntry pentry = {.seq = seq+cbytes, .pindex = 0,
                        .orient = FORWARD, .plen = plen};
    binary_seq_from_str(paths[i], plen, pentry.seq);
    pentrybuf_add(&set->members, pentry);
    set->seqs.len += cbytes + (plen+3)/4;
  }
}

static void _set_paths_colsets(SortedPathSet *set, const uint8_t *colset)
{
  size_t i;
  uint8_t *cset;
  for(i = 0; i < set->members.len; i++)
  {
    cset = sorted_path_colset(&set->members.data[i], set);
    memcpy(cset, colset, set->cbytes);
    colset += set->cbytes;
  }
}

static void _check_paths_colsets(const SortedPathSet *set, const uint8_t *colset)
{
  size_t i;
  const uint8_t *cset;
  for(i = 0; i < set->members.len; i++)
  {
    cset = sorted_path_colset(&set->members.data[i], set);
    TASSERT(memcmp(cset, colset, set->cbytes) == 0);
    colset += set->cbytes;
  }
}

// Construct a fake sorted_path_set and check that it paths are sorted
static void _check_fake_test_case(SortedPathSet *set, StrBuf *sbuf,
                                 const char **paths, size_t npaths)
{
  const PathEntry *entry;
  size_t i;

  TASSERT2(set->members.len == npaths, "%zu vs %zu", set->members.len, npaths);

  for(i = 0; i < set->members.len; i++) {
    entry = &set->members.data[i];
    strbuf_ensure_capacity(sbuf, entry->plen);
    binary_seq_to_str(entry->seq, entry->plen, sbuf->buff);
    TASSERT2(!strcmp(sbuf->buff,paths[i]), "'%s' vs '%s'", sbuf->buff, paths[i]);
  }
}

static void _fake_set_test()
{
  test_status("Testing sorted_path_set.c with fake sets...");

  StrBuf sbuf;
  strbuf_alloc(&sbuf, 1024);

  SortedPathSet set0, set1, set2;
  sorted_path_set_alloc(&set0);
  sorted_path_set_alloc(&set1);
  sorted_path_set_alloc(&set2);

  // Check sorted paths are sorted correctly

  const char *paths0a[4] = {"CG","C","T","CGC"}; // unsorted input
  const char *paths0b[4] = {"C","CG","CGC","T"};
  const char *paths0c[2] = {"CGC","T"};
  const char *paths0d[1] = {"CGC"};
  _fake_sorted_path_set(&set0, 1, paths0a, 4);
  _check_fake_test_case(&set0, &sbuf, paths0b, 4);
  sorted_path_set_slim(&set0);
  _check_fake_test_case(&set0, &sbuf, paths0c, 2);

  const char *paths1a[3] = {"TT","TT","A"}; // unsorted input
  const char *paths1b[3] = {"A","TT","TT"};
  const char *paths1c[2] = {"A","TT"};
  _fake_sorted_path_set(&set1, 1, paths1a, 3);
  _check_fake_test_case(&set1, &sbuf, paths1b, 3);

  // Filter set1 {A,TT,TT} against set0 {CGC,T} (no change)
  sorted_path_set_merge(&set0, &set1, false, NULL);
  _check_fake_test_case(&set1, &sbuf, paths1b, 3);

  // Filter set0 {CGC,T} against set1 {A,TT,TT}, including substring matches
  // removes T from set0 because it is a substr of TT
  sorted_path_set_merge(&set1, &set0, true, NULL);
  _check_fake_test_case(&set0, &sbuf, paths0d, 1);

  sorted_path_set_slim(&set1);
  _check_fake_test_case(&set1, &sbuf, paths1c, 2);

  const char *paths2a[1] = {"A"};
  const char *paths2b[1] = {"A"}; // just use as empty set
  _fake_sorted_path_set(&set2, 1, paths2a, 0);
  _check_fake_test_case(&set2, &sbuf, paths2b, 0);
  sorted_path_set_slim(&set2);
  _check_fake_test_case(&set2, &sbuf, paths2b, 0);

  sorted_path_set_dealloc(&set0);
  sorted_path_set_dealloc(&set1);
  sorted_path_set_dealloc(&set2);

  //
  // Slimming a set with colourset
  //

  // Set: [colourset,path]
  //   5=101 A
  //   6=110 AC
  //   5=101 ACG
  //   7=111 G
  //   0=000 T
  //   4=100 TA
  // Should go to:
  //   2=010 AC
  //   5=101 ACG
  //   7=111 G
  //   4=100 TA

  SortedPathSet set3, set4;
  sorted_path_set_alloc(&set3);
  sorted_path_set_alloc(&set4);

  const char *paths3a[6] = {"A","AC","ACG","G","T","TA"}; // unsorted input
  const uint8_t colset3a[6] = {5,6,5,7,0,4};

  const char *paths3b[6] = {"AC","ACG","G","TA"}; // output
  const uint8_t colset3b[6] = {2,5,7,4};

  // Create two copies of test set
  _fake_sorted_path_set(&set3, 1, paths3a, 6);
  _set_paths_colsets(&set3, colset3a);
  _fake_sorted_path_set(&set4, 1, paths3a, 6);
  _set_paths_colsets(&set4, colset3a);

  // Test slimming set3
  sorted_path_set_slim(&set3);
  _check_fake_test_case(&set3, &sbuf, paths3b, 4);
  _check_paths_colsets(&set3, colset3b);

  // Test slimming set3 x 2
  sorted_path_set_slim(&set3);
  _check_fake_test_case(&set3, &sbuf, paths3b, 4);
  _check_paths_colsets(&set3, colset3b);

  // Now attempt to merge set4 into the slimmed set3
  sorted_path_set_merge(&set3, &set4, true, NULL);
  // No change to set3
  _check_fake_test_case(&set3, &sbuf, paths3b, 4);
  _check_paths_colsets(&set3, colset3b);
  // Empty set4
  _check_fake_test_case(&set4, &sbuf, paths3b, 0);

  // Recreate set4
  _fake_sorted_path_set(&set4, 1, paths3a, 6);
  _set_paths_colsets(&set4, colset3a);

  // Set4:
  //   5=101 A
  //   6=110 AC
  //   5=101 ACG
  //   7=111 G
  //   0=000 T
  //   4=100 TA
  // Set3:
  //   2=010 AC
  //   5=101 ACG
  //   7=111 G
  //   4=100 TA

  // Merge set3 into set4
  sorted_path_set_merge(&set4, &set3, true, NULL);

  // Should have no effect to set4. should empty set4
  _check_fake_test_case(&set3, &sbuf, paths3b, 0);
  _set_paths_colsets(&set3, colset3b);
  _check_fake_test_case(&set4, &sbuf, paths3a, 6);
  _set_paths_colsets(&set4, colset3a);

  sorted_path_set_dealloc(&set3);
  sorted_path_set_dealloc(&set4);

  strbuf_dealloc(&sbuf);
}

//
// Real Graph
//

// Check we get paths in the order given
static void _test_real_sorted_path_set(dBGraph *graph, BinaryKmer bkmer,
                                       const char **paths, size_t npaths)
{
  SortedPathSet set;
  size_t i;
  dBNode node = db_graph_find(graph, bkmer);
  StrBuf sbuf;
  strbuf_alloc(&sbuf, 1024);

  sorted_path_set_alloc(&set);
  sorted_path_set_init(&set, &graph->pstore, node.key);

  TASSERT2(set.members.len == npaths, "%zu vs %zu", set.members.len, npaths);

  for(i = 0; i < npaths; i++)
  {
    PathEntry *entry = &set.members.data[i];
    strbuf_ensure_capacity(&sbuf, entry->plen);
    binary_seq_to_str(entry->seq, entry->plen, sbuf.buff);
    TASSERT2(strcmp(sbuf.buff, paths[i]) == 0, " %zu: '%s' vs '%s'",
             i, sbuf.buff, paths[i]);
  }

  strbuf_dealloc(&sbuf);
  sorted_path_set_dealloc(&set);
}

static void _real_graph_test()
{
  test_status("Testing sorted_path_set.c with real graph...");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1, path_max_mem = 1024;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 1024);
  // Graph data
  graph.bktlocks = calloc2(roundup_bits2bytes(graph.ht.num_of_buckets), 1);
  graph.col_edges = calloc2(graph.ht.capacity * ncols, sizeof(Edges));
  graph.col_covgs = calloc2(graph.ht.capacity * ncols, sizeof(Covg));
  graph.node_in_cols = calloc2(roundup_bits2bytes(graph.ht.capacity) * ncols, 1);

  // Path data
  path_store_alloc(&graph.pstore, path_max_mem, true, graph.ht.capacity, ncols);
  graph.pstore.kmer_locks = calloc2(roundup_bits2bytes(graph.ht.capacity), 1);

  // junctions:  >     >           <     <     <
  char seq0[] = "CCTGGGTGCGAATGACACCAAATCGAATGAC"; // a->d
  char seq1[] = "ACTGGGTGCGAATGACACCAAATCGAATGAT"; // b->e
  char seq2[] = "GACTATAGCGAATGACACCAAATCAGGGAGA"; // c->f
  char seq3[] = "GACTATAGCGAATGACACTTCTACCTGTCTC"; // c->g

  char s0r0[] = "CCTGGGTGCGAATGACACC"; // a->d
  char s0r1[] = "CCTGGGTGCGAATGACACCAAATCGAA"; // a->d
  char s0r2[] = "CCTGGGTGCGAATGACACCAAATCGAATGAC"; // a->d

  char s3r0[] = "CCTGGGTGCGAATGACACT"; // a->g

  // a--+--+--+--+--+--d
  // b-/  /    \  \  \_e
  // c___/      \  \___f
  //             \_____g

  build_graph_from_str_mt(&graph, 0, seq0, strlen(seq0));
  build_graph_from_str_mt(&graph, 0, seq1, strlen(seq1));
  build_graph_from_str_mt(&graph, 0, seq2, strlen(seq2));
  build_graph_from_str_mt(&graph, 0, seq3, strlen(seq3));

  // Set up alignment correction params
  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .ins_gap_min = 0, .ins_gap_max = 0,
                            .one_way_gap_traverse = true, .use_end_check = true,
                            .max_context = 10,
                            .gap_variance = 0.1, .gap_wiggle = 5};

  AsyncIOReadTask io = {.file1 = NULL, .file2 = NULL,
                        .fq_offset = 0, .interleaved = false};

  // Load paths
  CorrectAlnReadsTask task = {.files = io, .fq_cutoff = 0, .hp_cutoff = 0,
                              .matedir = READPAIR_FR, .crt_params = params,
                              .ptr = NULL};

  AsyncIOData iodata;
  asynciodata_alloc(&iodata);

  size_t nworkers = 1;
  GenPathWorker *wrkrs = gen_paths_workers_alloc(nworkers, &graph, NULL);

  // Add paths
  _test_add_paths(&graph, &iodata, &task, wrkrs, s3r0, 3, 3, 3); // T
  _test_add_paths(&graph, &iodata, &task, wrkrs, s0r0, 3, 1, 3); // C
  _test_add_paths(&graph, &iodata, &task, wrkrs, s0r2, 4, 2, 4); // CGC
  _test_add_paths(&graph, &iodata, &task, wrkrs, s0r1, 2, 0, 2); // CG

  // CCTGGGTGCGA CCTGGGTGCGA:0 C,CG,CGC,T [a->d][a->g]
  // TGCGAATGACA TGCGAATGACA:0 C,CG,CGC,T [a->d][a->g]
  // GGTGTCATTCG CGAATGACACC:1 AG [d->a]
  // CGATTTGGTGT ACACCAAATCG:1 AG [d->a]
  // GTCATTCGATT AATCGAATGAC:1 AG [d->a]
  // AGTGTCATTCG AGTGTCATTCG:0 AG [g->a]

  // Test path store
  _test_path_store(&graph);

  // Test sorted_path_set
  BinaryKmer bkmer;
  const char *seqs0[4] = {"C","CG","CGC", "T"};
  const char *seqs1[1] = {"AG"};

  bkmer = binary_kmer_from_str("CCTGGGTGCGA", kmer_size);
  _test_real_sorted_path_set(&graph, bkmer, seqs0, 4);
  bkmer = binary_kmer_from_str("TGCGAATGACA", kmer_size);
  _test_real_sorted_path_set(&graph, bkmer, seqs0, 4);
  bkmer = binary_kmer_from_str("GGTGTCATTCG", kmer_size);
  _test_real_sorted_path_set(&graph, bkmer, seqs1, 1);
  bkmer = binary_kmer_from_str("CGATTTGGTGT", kmer_size);
  _test_real_sorted_path_set(&graph, bkmer, seqs1, 1);
  bkmer = binary_kmer_from_str("GTCATTCGATT", kmer_size);
  _test_real_sorted_path_set(&graph, bkmer, seqs1, 1);
  bkmer = binary_kmer_from_str("AGTGTCATTCG", kmer_size);
  _test_real_sorted_path_set(&graph, bkmer, seqs1, 1);

  // Test node with no paths
  bkmer = binary_kmer_from_str("GCGAATGACAC", kmer_size);
  _test_real_sorted_path_set(&graph, bkmer, NULL, 0);

  gen_paths_workers_dealloc(wrkrs, nworkers);

  ctx_free(graph.bktlocks);
  ctx_free(graph.node_in_cols);
  ctx_free(graph.col_edges);
  ctx_free(graph.col_covgs);

  asynciodata_dealloc(&iodata);
  path_store_dealloc(&graph.pstore);
  db_graph_dealloc(&graph);
}

void test_sorted_paths()
{
  _real_graph_test();
  _fake_set_test();
}
