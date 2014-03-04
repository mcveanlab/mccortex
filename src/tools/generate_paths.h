#ifndef GENERATE_PATHS_H_
#define GENERATE_PATHS_H_

#include "seq_file.h"

#include "cortex_types.h"
#include "db_graph.h"
#include "loading_stats.h"
#include "correct_reads_input.h"

typedef struct GenPathWorker GenPathWorker;

// Estimate memory required per worker thread
size_t gen_paths_worker_est_mem(const dBGraph *db_graph);

GenPathWorker* gen_paths_workers_alloc(size_t n, dBGraph *graph);
void gen_paths_workers_dealloc(GenPathWorker *mem, size_t n);

// Add a single contig using a given worker
void gen_path_worker_seq(GenPathWorker *wrkr, AsyncIOData *data,
                         const CorrectAlnReadsTask *task);

// workers array must be at least as long as tasks
void generate_paths(CorrectAlnReadsTask *tasks, size_t num_tasks,
                    GenPathWorker *workers, size_t num_workers);

// Save gap size distribution
// base_fmt is the beginning of the file name - the reset is <num>.csv or something
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void gen_paths_dump_gap_sizes(const char *path,
                              const uint64_t *arr, size_t arrlen,
                              size_t kmer_size, bool insert_sizes,
                              size_t nreads);

// Get histogram array
const uint64_t* gen_paths_get_ins_gap(GenPathWorker *worker, size_t *len);
const uint64_t* gen_paths_get_err_gap(GenPathWorker *worker, size_t *len);

void gen_paths_get_stats(const GenPathWorker *worker, size_t num_workers,
                         LoadingStats *stats);

#endif /* GENERATE_PATHS_H_ */
