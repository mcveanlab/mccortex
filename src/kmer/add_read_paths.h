#ifndef ADD_READ_PATH_H_
#define ADD_READ_PATH_H_

#include "file_reader.h"

typedef struct PathsWorkerPool PathsWorkerPool;

PathsWorkerPool* paths_worker_pool_new(size_t num_of_threads,
                                       dBGraph *db_graph, size_t max_gap_limit);

void paths_worker_pool_dealloc(PathsWorkerPool *pool);

void add_read_paths_to_graph(PathsWorkerPool *pool,
                             seq_file_t *sf1, seq_file_t *sf2,
                             size_t gap_limit, size_t ctx_col, size_t ctp_col,
                             const SeqLoadingPrefs *prefs,
                             SeqLoadingStats *stats);

#endif /* ADD_READ_PATH_H_ */
