#include "global.h"
#include "all_tests.h"
#include "db_graph.h"

// Common functions here
FILE *ctx_tst_out = NULL;

size_t tests_num_run = 0, tests_num_failed = 0;

//
// Useful functions
//
void fill_rand_bytes(uint8_t *arr, size_t n)
{
  uint32_t r; size_t i;
  for(i = 0; i+4 <= n; i+=4) { r = rand(); memcpy(&arr[i], &r, 4); }
  r = rand();
  for(; i < n; i++) { arr[i] = (uint8_t)(r&0xff); r >>= 8; }
}

void rand_nucs(Nucleotide *nucs, size_t len)
{
  if(!len) return;
  size_t i, r = 0;
  for(i = 0; i < len; i++) {
    if((i & 15) == 0) r = (size_t)rand(); // 2 bits per cycle, 32 bits in rand()
    nucs[i] = r&3;
    r >>= 2;
  }
}

void bitarr_tostr(const uint8_t *arr, size_t len, char *str)
{
  size_t i, j;
  for(i = len-1; i != SIZE_MAX; i--) {
    for(j = 7; j != SIZE_MAX; j--) {
      *str = (arr[i]&(1U<<j)) ? '1' : '0';
      str++;
    }
  }
  *str = '\0';
}

//
// Graph setup
//

void _construct_graph_with_paths(dBGraph *graph,
                                 size_t kmer_size, size_t ncols,
                                 char **seqs, size_t nseqs,
                                 CorrectAlnParam path_params)
{
  size_t i;
  db_graph_alloc(graph, kmer_size, ncols, ncols, 1024);

  // Graph data
  graph->bktlocks = calloc2(roundup_bits2bytes(graph->ht.num_of_buckets), 1);
  graph->col_edges = calloc2(graph->ht.capacity * ncols, sizeof(Edges));
  graph->col_covgs = calloc2(graph->ht.capacity * ncols, sizeof(Covg));
  graph->node_in_cols = calloc2(roundup_bits2bytes(graph->ht.capacity) * ncols, 1);

  // Path data
  path_store_alloc(&graph->pstore, 1024, true, graph->ht.capacity, ncols);
  graph->pstore.kmer_locks = calloc2(roundup_bits2bytes(graph->ht.capacity), 1);

  // Build graph
  for(i = 0; i < nseqs; i++)
    build_graph_from_str_mt(graph, 0, seqs[i], strlen(seqs[i]));

  GenPathWorker *gen_path_wrkr = gen_paths_workers_alloc(1, graph, NULL);

  for(i = 0; i < nseqs; i++)
    gen_paths_from_str_mt(gen_path_wrkr, seqs[i], path_params);

  gen_paths_workers_dealloc(gen_path_wrkr, 1);
}

void _deconstruct_graph_with_paths(dBGraph *graph)
{
  free(graph->bktlocks);
  free(graph->node_in_cols);
  free(graph->col_edges);
  free(graph->col_covgs);

  path_store_dealloc(&graph->pstore);
  db_graph_dealloc(graph);
}
