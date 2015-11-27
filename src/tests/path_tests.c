#include "global.h"
#include "all_tests.h"
#include "db_graph.h"
#include "build_graph.h"
#include "generate_paths.h"
#include "gpath_checks.h"

//       junctions:  >     >           <     <     <
const char seq0[] = "CCTGGGTGCGAATGACACCAAATCGAATGAC"; // a->d
const char seq1[] = "ACTGGGTGCGAATGACACCAAATCGAATGAT"; // b->e
const char seq2[] = "GACTATAGCGAATGACACCAAATCAGGGAGA"; // c->f
const char seq3[] = "GACTATAGCGAATGACACTTCTACCTGTCTC"; // c->g

// a--+--+--+--+--+--d
// b-/  /    \  \  \_e
// c___/      \  \___f
//             \_____g

// "CCTGGGTGCGA", "CCTGGGTGCGAATGACACCAAATCGAATGAC"
// "ACTGGGTGCGA", "ACTGGGTGCGAATGACACCAAATCGAATGAT"

const char kmerA[] = "CCTGGGTGCGA"; // a->
#define NPATHS_A 1
const char *kmerApaths[NPATHS_A] = {"CCTGGGTGCGAATGACACCAAATCGAATGAC"};

const char kmerB[] = "ACTGGGTGCGA"; // b->
#define NPATHS_B 1
const char *kmerBpaths[NPATHS_B] = {"ACTGGGTGCGAATGACACCAAATCGAATGAT"};

// where seq0 and seq1 meet
const char kmerAB[] = "TGCGAATGACA"; // {a->,b->}->
#define NPATHS_AB 2
const char *kmerABpaths[NPATHS_AB] = {"TGCGAATGACACCAAATCGAATGAC",
                                      "TGCGAATGACACCAAATCGAATGAT"};

const char kmerC[] = "AGCGAATGACA"; // c->
#define NPATHS_C 2
const char *kmerCpaths[NPATHS_C] = {"AGCGAATGACACCAAATCA", // c->f
                                     "AGCGAATGACACT"};      // c->g

/* Reverse */
const char kmerG[] = "AGTGTCATTCG"; // ->g revcmp(CGAATGACACT)
#define NPATHS_G 1
const char *kmerGpaths[NPATHS_G] = {"AGTGTCATTCGCT"}; // revcmp(AGCGAATGACACT)

const char kmerF[] = "TGATTTGGTGT"; // ->f revcmp(ACACCAAATCA)
#define NPATHS_F 1
const char *kmerFpaths[NPATHS_F] = {"TGATTTGGTGTCATTCGCT"}; // revcmp(AGCGAATGACACCAAATCA)

const char kmerE[] = "ATCATTCGATT"; // ->e revcmp(AATCGAATGAT)
#define NPATHS_E 1
const char *kmerEpaths[NPATHS_E] = {"ATCATTCGATTTGGTGTCATTCGCACCCAGT"}; // revcmp(ACTGGGTGCGAATGACACCAAATCGAATGAT)

const char kmerD[] = "GTCATTCGATT"; // ->d revcmp(AATCGAATGAC)
#define NPATHS_D 1
const char *kmerDpaths[NPATHS_D] = {"GTCATTCGATTTGGTGTCATTCGCACCCAGG"}; // revcmp(CCTGGGTGCGAATGACACCAAATCGAATGAC)

const char kmerDEF[] = "GGTGTCATTCG"; // ->def revcmp(CGAATGACACC)
#define NPATHS_DEF 3
const char *kmerDEFpaths[NPATHS_DEF] = {"GGTGTCATTCGCACCCAGG", // revcmp(CCTGGGTGCGAATGACACC)
                                        "GGTGTCATTCGCACCCAGT", // revcmp(ACTGGGTGCGAATGACACC)
                                        "GGTGTCATTCGCT"};      // revcmp(AGCGAATGACACC)

const char kmerDE[]  = "CGATTTGGTGT"; // ->de  revcmp(ACACCAAATCG)
#define NPATHS_DE 2
const char *kmerDEpaths[NPATHS_DE] = {"CGATTTGGTGTCATTCGCACCCAGG", // revcmp(CCTGGGTGCGAATGACACCAAATCG)
                                      "CGATTTGGTGTCATTCGCACCCAGT"};// revcmp(ACTGGGTGCGAATGACACCAAATCG)

static void _check_node_paths(const char *kmer,
                              const char **path_strs, size_t npaths,
                              size_t colour, const dBGraph *graph)
{
  TASSERT(strlen(kmer) == graph->kmer_size);

  const GPath *paths[npaths]; // corresponding to path_strs
  memset(paths, 0, sizeof(paths));
  size_t i, num_paths_seen = 0;

  const GPathStore *gpstore = &graph->gpstore;
  dBNode node = db_graph_find_str(graph, kmer);

  const GPath *path = gpath_store_fetch_traverse(gpstore, node.key);
  dBNodeBuffer nbuf;
  SizeBuffer jposbuf;
  db_node_buf_alloc(&nbuf, 64);
  size_buf_alloc(&jposbuf, 64);

  #define MAX_SEQ 128
  char seq[MAX_SEQ];

  for(; path != NULL; path = path->next)
  {
    if(path->orient == node.orient &&
       gpath_has_colour(path, gpstore->gpset.ncols, colour))
    {
      TASSERT(num_paths_seen < npaths);
      db_node_buf_reset(&nbuf);
      gpath_fetch(node, path, &nbuf, &jposbuf, colour, graph);
      if(nbuf.len > MAX_SEQ) die("Too many nodes. Cannot continue. %zu", nbuf.len);
      db_nodes_to_str(nbuf.b, nbuf.len, graph, seq);
      TASSERT(strlen(seq) == graph->kmer_size + nbuf.len - 1);
      for(i = 0; i < npaths; i++) {
        if(strcmp(path_strs[i],seq) == 0) {
          TASSERT(paths[i] == NULL, "Duplicate paths: %s", seq);
          paths[i] = path;
          break;
        }
      }
      TASSERT2(i < npaths, "Path not found: %s", seq);
      num_paths_seen++;
    }
  }

  TASSERT(num_paths_seen == npaths);

  for(i = 0; i < npaths; i++) {
    TASSERT2(paths[i] != NULL, "path not in graph: %s", path_strs[i]);
  }

  db_node_buf_dealloc(&nbuf);
  size_buf_dealloc(&jposbuf);
}

static void _test_add_paths()
{
  test_status("Testing adding paths in generate_paths.c and gpath_fetch()");

  // Construct 1 colour graph with kmer-size=11
  dBGraph graph;
  size_t kmer_size = 11, ncols = 1;

  db_graph_alloc(&graph, kmer_size, ncols, ncols, 1024,
                 DBG_ALLOC_EDGES | DBG_ALLOC_COVGS |
                 DBG_ALLOC_BKTLOCKS | DBG_ALLOC_NODE_IN_COL);

  // Create a path store that tracks path counts
  gpath_store_alloc(&graph.gpstore,
                    graph.num_of_cols, graph.ht.capacity,
                    0, ONE_MEGABYTE, true, false);

  // Create path hash table for fast lookup
  gpath_hash_alloc(&graph.gphash, &graph.gpstore, ONE_MEGABYTE);

  build_graph_from_str_mt(&graph, 0, seq0, strlen(seq0), false);
  build_graph_from_str_mt(&graph, 0, seq1, strlen(seq1), false);
  build_graph_from_str_mt(&graph, 0, seq2, strlen(seq2), false);
  build_graph_from_str_mt(&graph, 0, seq3, strlen(seq3), false);

  // Set up alignment correction params
  CorrectAlnParam params = {.ctpcol = 0, .ctxcol = 0,
                            .frag_len_min = 0, .frag_len_max = 0,
                            .one_way_gap_traverse = true, .use_end_check = true,
                            .max_context = 10,
                            .gap_variance = 0.1, .gap_wiggle = 5};

  all_tests_add_paths(&graph, seq0, params, 5, 5); // path lens: 3+3+2+2+2
  all_tests_add_paths(&graph, seq1, params, 5, 2); // path lens: 3+3+2+2+2
  all_tests_add_paths(&graph, seq2, params, 3, 2); // path lens: 1+1+1
  all_tests_add_paths(&graph, seq3, params, 2, 1); // path lens: 1+1

  // Test path store
  gpath_checks_all_paths(&graph, 1); // use one thread

  // Test path content
  _check_node_paths(kmerA,  kmerApaths,  NPATHS_A,  0, &graph);
  _check_node_paths(kmerB,  kmerBpaths,  NPATHS_B,  0, &graph);
  _check_node_paths(kmerAB, kmerABpaths, NPATHS_AB, 0, &graph);
  _check_node_paths(kmerC,  kmerCpaths,  NPATHS_C,  0, &graph);
  _check_node_paths(kmerG,  kmerGpaths,  NPATHS_G,  0, &graph);
  _check_node_paths(kmerF,  kmerFpaths,  NPATHS_F,  0, &graph);
  _check_node_paths(kmerE,  kmerEpaths,  NPATHS_E,  0, &graph);
  _check_node_paths(kmerD,  kmerDpaths,  NPATHS_D,  0, &graph);
  _check_node_paths(kmerDEF,kmerDEFpaths,NPATHS_DEF,0, &graph);
  _check_node_paths(kmerDE, kmerDEpaths, NPATHS_DE, 0, &graph);

  db_graph_dealloc(&graph);
}

void test_paths()
{
  _test_add_paths();
}
