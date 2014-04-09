#include "global.h"
#include "all_tests.h"
#include "db_graph.h"
#include "build_graph.h"
#include "path_store.h"
#include "generate_paths.h"
#include "graph_paths.h"

void test_paths()
{
  test_status("Testing adding paths in generate_paths.c");

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

  _test_add_paths(&graph, &iodata, &task, wrkrs, seq0, 5, 5, 5); // path lens: 3+3+2+2+2
  _test_add_paths(&graph, &iodata, &task, wrkrs, seq1, 5, 2, 5); // path lens: 3+3+2+2+2
  _test_add_paths(&graph, &iodata, &task, wrkrs, seq2, 3, 2, 3); // path lens: 1+1+1
  _test_add_paths(&graph, &iodata, &task, wrkrs, seq3, 2, 1, 2); // path lens: 1+1

  path_store_combine_updated_paths(&graph.pstore);

  // DEV: Actually test path content, colours set etc
  // seq0 fw
  // CCTGGGTGCGA:0 len:3 col:0  CGC
  // TGCGAATGACA:0 len:3 col:0  CGC
  // seq0 rv
  // AATCGAATGAC:1 len:2 col:0  AG
  // ACACCAAATCG:1 len:2 col:0  AG
  // CGAATGACACC:1 len:2 col:0  AG

  // Test path store
  _test_path_store(&graph);

  gen_paths_workers_dealloc(wrkrs, nworkers);

  ctx_free(graph.bktlocks);
  ctx_free(graph.node_in_cols);
  ctx_free(graph.col_edges);
  ctx_free(graph.col_covgs);

  asynciodata_dealloc(&iodata);
  path_store_dealloc(&graph.pstore);
  db_graph_dealloc(&graph);
}

/*
typedef struct {
  size_t npaths, npath_cap, nbases, nbases_cap, *order;
  PathLen *len_orients;
  Nucleotide *bases;
} PathList;

static inline void path_list_alloc(PathList *plist)
{
  plist->npaths = plist->nbases = 0;
  plist->npath_cap = 512;
  plist->nbases_cap = 1024;
  plist->len_orients = malloc2(plist->npath_cap * sizeof(*plist->len_orients));
  plist->order = malloc2(plist->npath_cap * sizeof(*plist->order));
  plist->bases = malloc2(plist->nbases_cap * sizeof(*plist->bases));
}

static inline void path_list_init(PathList *plist) {
  plist->npaths = plist->nbases = 0;
}

static inline void path_list_dealloc(PathList *plist) {
  ctx_free(plist->order); ctx_free(plist->len_orients); ctx_free(plist->bases);
}

static inline void add_path(PathList *plist,
                            const PathStore *pstore, PathIndex pi)
{
  PathLen len, merged; Orientation orient;
  merged = packedpath_get_len_orient(pstore->store + pi, pstore->colset_bytes,
                                     &len, &orient);

  if(plist->nbases + len > plist->nbases_cap) {
    plist->nbases_cap = roundup2pow(plist->nbases + len);
    plist->bases = realloc2(plist->bases, plist->nbases_cap*sizeof(*plist->bases));
    plist->order = realloc2(plist->order, plist->nbases_cap*sizeof(*plist->order));
  }
  if(plist->npaths + 1 > plist->npath_cap) {
    plist->npath_cap *= 2;
    plist->len_orients = realloc2(plist->len_orients, plist->npath_cap);
  }

  plist->len_orients[plist->npaths++] = merged;
  // path_store_fetch_bases(pstore, pi, plist->bases+plist->nbases, len);
  plist->nbases += len;
}

static inline int plist_cmp(const void *a, const void *b, void *arg)
{
  size_t x = *(const size_t*)a, y = *(const size_t*)b;
  PathList *plist = (PathList*)arg;
  ptrdiff_t cmp = plist->len_orients[y] - plist->len_orients[x];
  if(cmp != 0) return cmp;
  // DEV:
  // return memcmp(plist->bases+);
  return 0;
}

static inline void compare_kmer_paths(hkey_t node,
                                      const dBGraph *dbg1, const dBGraph *dbg2,
                                      PathList *plist1,
                                      PathList *plist2)
{
  TASSERT(binary_kmers_are_equal(dbg1->ht.table[node], dbg2->ht.table[node]));

  // Fecth path list, sort, and compare paths
  if((db_node_paths(dbg1, node) == PATH_NULL) !=
     (db_node_paths(dbg2, node) == PATH_NULL))
  {
    char bstr[MAX_KMER_SIZE+1];
    binary_kmer_to_str(db_node_get_bkmer(dbg1, node), dbg1->kmer_size, bstr);
    die("Kmer has path in only one graph [%zu vs %zu]: %s",
        (size_t)db_node_paths(dbg1, node), (size_t)db_node_paths(dbg2, node), bstr);
  }

  path_list_init(plist1);
  path_list_init(plist2);

  PathIndex pi = db_node_paths(dbg1, node);
  while(pi != PATH_NULL) {
    add_path(plist1, &dbg1->pstore, pi);
    pi = packedpath_get_prev(dbg1->pstore.store+pi);
  }
  pi = db_node_paths(dbg2, node);
  while(pi != PATH_NULL) {
    add_path(plist2, &dbg2->pstore, pi);
    pi = packedpath_get_prev(dbg2->pstore.store+pi);
  }

  if(plist1->npaths != plist2->npaths) die("Mismatch in lengths");

  // Sort plist1, plist2
  size_t i;
  for(i = 0; i < plist1->npaths; i++) plist1->order[i] = plist2->order[i] = i;
  // sort_r(plist1->order, plist1->npaths, sizeof(*plist1->order), plist_cmp);

  // Compare
  // reset plist1, plist2
}

//
void test_path_stores_match(const dBGraph *dbg1, const dBGraph *dbg2)
{
  // Test not valid unless graphs were built in the same way
  TASSERT(dbg1->ht.capacity == dbg2->ht.capacity);
  TASSERT(dbg1->ht.num_of_buckets == dbg2->ht.num_of_buckets);
  TASSERT(dbg1->ht.bucket_size == dbg2->ht.bucket_size);
  TASSERT(dbg1->kmer_size == dbg2->kmer_size);

  PathList plist1, plist2;

  path_list_alloc(&plist1);
  path_list_alloc(&plist2);

  HASH_ITERATE(&dbg1->ht, compare_kmer_paths, dbg1, dbg2, &plist1, &plist2);

  path_list_dealloc(&plist1);
  path_list_dealloc(&plist2);
}
*/
