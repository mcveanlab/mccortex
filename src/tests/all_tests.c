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

void all_tests_add_paths_multi(dBGraph *graph, const char **seqs, size_t nseqs,
                               CorrectAlnParam params,
                               int exp_npaths, int exp_nkmers)
{
  size_t npaths = graph->gpstore.num_paths;
  size_t nkmers = graph->gpstore.num_kmers_with_paths;

  size_t i, nworkers = 1;
  GenPathWorker *wrkrs = gen_paths_workers_alloc(nworkers, graph);

  // Set up asyncio input data
  AsyncIOInput io = {.file1 = NULL, .file2 = NULL,
                     .fq_offset = 0, .interleaved = false};

  CorrectAlnInput task = {.files = io, .fq_cutoff = 0, .hp_cutoff = 0,
                          .matedir = READPAIR_FR, .crt_params = params,
                          .out_base = NULL, .output = NULL};

  AsyncIOData iodata;
  asynciodata_alloc(&iodata);
  seq_read_reset(&iodata.r2);
  iodata.fq_offset1 = iodata.fq_offset2 = 0;
  iodata.ptr = NULL;

  // Add paths
  for(i = 0; i < nseqs; i++) {
    seq_read_set(&iodata.r1, seqs[i]);
    gen_paths_worker_seq(wrkrs, &iodata, &task);
  }

  asynciodata_dealloc(&iodata);
  gen_paths_workers_dealloc(wrkrs, nworkers);

  // Check we added the right number of paths
  if(exp_npaths >= 0) {
    TASSERT2(graph->gpstore.num_paths == npaths + (size_t)exp_npaths, "%zu %zu %zu",
             (size_t)graph->gpstore.num_paths, (size_t)npaths, (size_t)exp_npaths);
  }

  if(exp_nkmers >= 0) {
    TASSERT(graph->gpstore.num_kmers_with_paths == nkmers + (size_t)exp_nkmers);
  }
}

void all_tests_add_paths(dBGraph *graph, const char *seq,
                         CorrectAlnParam params,
                         int exp_npaths, int exp_nkmers)
{
  all_tests_add_paths_multi(graph, &seq, 1, params, exp_npaths, exp_nkmers);
}

void all_tests_construct_graph(dBGraph *graph,
                               size_t kmer_size, size_t ncols,
                               const char **seqs, size_t nseqs,
                               CorrectAlnParam path_params)
{
  size_t i;
  db_graph_alloc(graph, kmer_size, ncols, ncols, 1024,
                 DBG_ALLOC_EDGES | DBG_ALLOC_COVGS |
                 DBG_ALLOC_NODE_IN_COL | DBG_ALLOC_BKTLOCKS);

  // Path data
  gpath_store_alloc(&graph->gpstore, ncols, graph->ht.capacity,
                    0, ONE_MEGABYTE, true, false);

  // Don't use links to add new links
  gpath_store_split_read_write(&graph->gpstore);

  // Allocate path hash table just in case
  gpath_hash_alloc(&graph->gphash, &graph->gpstore, ONE_MEGABYTE);

  // Build graph
  for(i = 0; i < nseqs; i++)
    build_graph_from_str_mt(graph, 0, seqs[i], strlen(seqs[i]), false);

  gpath_store_merge_read_write(&graph->gpstore);

  graph->num_of_cols_used = MAX2(graph->num_of_cols_used, 1);

  all_tests_add_paths_multi(graph, seqs, nseqs, path_params, -1, -1);
}
