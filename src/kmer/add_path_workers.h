#ifndef ADD_PATH_WORKERS_H_
#define ADD_PATH_WORKERS_H_

#include "file_reader.h"

typedef struct PathsWorkerPool PathsWorkerPool;

PathsWorkerPool* path_workers_pool_new(size_t num_of_threads,
                                       dBGraph *db_graph,
                                       size_t max_gap_limit);

// Read numbers are for printing out only
void path_workers_pool_dealloc(PathsWorkerPool *pool,
                               size_t num_se_reads, size_t num_pe_readpairs);

void path_workers_wait_til_finished(PathsWorkerPool *pool);

void path_workers_add_paths_to_graph(PathsWorkerPool *pool,
                                     seq_file_t *sf1, seq_file_t *sf2,
                                     size_t ctx_col, size_t ctp_col,
                                     size_t ins_gap_min, size_t ins_gap_max,
                                     const SeqLoadingPrefs *prefs,
                                     SeqLoadingStats *stats);

#endif
