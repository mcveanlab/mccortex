#include "global.h"
#include "all_tests.h"
#include "db_graph.h"
#include "dna.h"

// Common functions here
FILE *ctx_tst_out = NULL;

size_t tests_num_run = 0, tests_num_failed = 0;

//
// Useful functions
//
void rand_bytes(uint8_t *arr, size_t n)
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

void rand_bases(char *bases, size_t len)
{
  if(!len) return;
  size_t i, r = 0;
  for(i = 0; i < len; i++) {
    if((i & 15) == 0) r = (size_t)rand(); // 2 bits per cycle, 32 bits in rand()
    bases[i] = dna_nuc_to_char(r&3);
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
  graph->bktlocks = ctx_calloc(roundup_bits2bytes(graph->ht.num_of_buckets), 1);
  graph->col_edges = ctx_calloc(graph->ht.capacity * ncols, sizeof(Edges));
  graph->col_covgs = ctx_calloc(graph->ht.capacity * ncols, sizeof(Covg));
  graph->node_in_cols = ctx_calloc(roundup_bits2bytes(graph->ht.capacity) * ncols, 1);

  // Path data
  gpath_store_alloc(&graph->gpstore, ncols, graph->ht.capacity,
                    ONE_MEGABYTE, true, false);

  // Allocate path hash table just in case
  gpath_hash_alloc(&graph->gphash, &graph->gpstore, ONE_MEGABYTE);

  // Build graph
  for(i = 0; i < nseqs; i++)
    build_graph_from_str_mt(graph, 0, seqs[i], strlen(seqs[i]));

  graph->num_of_cols_used = MAX2(graph->num_of_cols_used, 1);

  GenPathWorker *gen_path_wrkr = gen_paths_workers_alloc(1, graph);

  for(i = 0; i < nseqs; i++)
    gen_paths_from_str_mt(gen_path_wrkr, seqs[i], path_params);

  gen_paths_workers_dealloc(gen_path_wrkr, 1);
}

void _test_add_paths(dBGraph *graph,
                     AsyncIOData *iodata, CorrectAlnInput *task,
                     GenPathWorker *wrkrs, char *seq,
                     size_t exp_npaths, size_t exp_nkmers)
{
  size_t npaths = graph->gpstore.num_paths;
  size_t nkmers = graph->gpstore.num_kmers_with_paths;

  // Set up asyncio input data
  seq_read_set(&iodata->r1, seq);
  seq_read_reset(&iodata->r2);
  iodata->fq_offset1 = iodata->fq_offset2 = 0;
  iodata->ptr = NULL;
  // Add paths
  gen_paths_worker_seq(wrkrs, iodata, task);

  // Check we added the right number of paths
  TASSERT2(graph->gpstore.num_paths == npaths + exp_npaths, "%zu %zu %zu",
           (size_t)graph->gpstore.num_paths, (size_t)npaths, (size_t)exp_npaths);
  TASSERT(graph->gpstore.num_kmers_with_paths == nkmers + exp_nkmers);
}
