#ifndef ADD_READ_PATH_H_
#define ADD_READ_PATH_H_

#include "file_reader.h"
#include "graph_walker.h"
#include "repeat_walker.h"

// typedef struct PathsWorkerPool PathsWorkerPool;

// PathsWorkerPool* paths_worker_pool_new(size_t num_of_threads,
//                                        dBGraph *db_graph, size_t max_gap_limit);

// void paths_worker_pool_dealloc(PathsWorkerPool *pool);

// void wait_for_jobs_to_finish(PathsWorkerPool *pool);

// void add_read_paths_to_graph(PathsWorkerPool *pool,
//                              seq_file_t *sf1, seq_file_t *sf2,
//                              size_t gap_limit, size_t ctx_col, size_t ctp_col,
//                              const SeqLoadingPrefs *prefs,
//                              SeqLoadingStats *stats);


typedef struct {
  read_t r1, r2;
  Colour ctp_col, ctx_col;
  size_t gap_limit;
  int qcutoff1, qcutoff2;
  int hp_cutoff;
} AddPathsJob;


void add_read_paths_init();
void add_read_paths_cleanup();

void add_read_paths(const AddPathsJob *job, dBNodeBuffer *nodebuf,
                    GraphWalker *wlk, RepeatWalker *rptwlk,
                    uint64_t *insert_sizes, uint64_t *gap_sizes,
                    dBGraph *db_graph);

#endif /* ADD_READ_PATH_H_ */
