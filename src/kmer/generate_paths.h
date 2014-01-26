#ifndef GENERATE_PATHS_H_
#define GENERATE_PATHS_H_

#include "seq_file.h"

#include "cortex_types.h"
#include "db_graph.h"
#include "loading_stats.h"

typedef struct
{
  seq_file_t *const file1, *const file2;
  const Colour ctpcol, ctxcol;
  const uint32_t ins_gap_min, ins_gap_max;
  const uint8_t fq_offset, fq_cutoff, hp_cutoff;
  const boolean read_pair_FR; // set to false if reads are FF
  const boolean one_way_gap_traverse; // set to false for more error prone algo
} GeneratePathsTask;

typedef struct GenPathWorker GenPathWorker;

extern boolean gen_paths_print_contigs;

// Estimate memory required per worker thread
size_t gen_paths_worker_est_mem(const dBGraph *db_graph);

GenPathWorker* gen_paths_workers_alloc(size_t n, dBGraph *graph);
void gen_paths_workers_dealloc(GenPathWorker *mem, size_t n);

// Add a single contig using a given worker
void gen_path_worker_seq(GenPathWorker *wrkr, const GeneratePathsTask *task,
                         const char *seq, size_t len);

// workers array must be at least as long as tasks
void generate_paths(GeneratePathsTask *tasks, size_t num_tasks,
                    GenPathWorker *workers, size_t num_workers);

// Save gap size distribution
// base_fmt is the beginning of the file name - the reset is <num>.csv or something
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void gen_paths_dump_gap_sizes(const char *base_fmt,
                              const uint64_t *arr, size_t arrlen,
                              size_t kmer_size, boolean insert_sizes,
                              size_t nreads);

// Get histogram array
const uint64_t* gen_paths_get_ins_gap(GenPathWorker *worker, size_t *len);
const uint64_t* gen_paths_get_err_gap(GenPathWorker *worker, size_t *len);

void gen_paths_get_stats(GenPathWorker *worker, size_t num_workers,
                         SeqLoadingStats *stats);

#endif /* GENERATE_PATHS_H_ */
