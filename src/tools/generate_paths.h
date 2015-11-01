#ifndef GENERATE_PATHS_H_
#define GENERATE_PATHS_H_

#include "seq_file.h"

#include "cortex_types.h"
#include "db_graph.h"
#include "seq_loading_stats.h"
#include "correct_aln_input.h"

typedef struct GenPathWorker GenPathWorker;

// Estimate memory required per worker thread
size_t gen_paths_worker_est_mem(const dBGraph *db_graph);

GenPathWorker* gen_paths_workers_alloc(size_t n, dBGraph *graph);

void gen_paths_workers_dealloc(GenPathWorker *mem, size_t n);

// Add a single contig using a given worker
void gen_paths_worker_seq(GenPathWorker *wrkr, AsyncIOData *data,
                          const CorrectAlnInput *task);

// For testing
void gen_paths_from_str_mt(GenPathWorker *gen_path_wrkr, char *seq,
                           CorrectAlnParam params);

// workers array must be at least as long as tasks
void generate_paths(CorrectAlnInput *tasks, size_t num_tasks,
                    GenPathWorker *workers, size_t num_workers);

CorrectAlnStats* gen_paths_get_aln_stats(GenPathWorker *wrkr);
SeqLoadingStats* gen_paths_get_stats(GenPathWorker *wrkr);

#endif /* GENERATE_PATHS_H_ */
