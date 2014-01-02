#ifndef ADD_READ_PATH_H_
#define ADD_READ_PATH_H_

#include "file_reader.h"
#include "graph_walker.h"
#include "repeat_walker.h"

typedef struct {
  read_t r1, r2;
  Colour ctp_col, ctx_col;
  size_t ins_gap_min, ins_gap_max;
  int qcutoff1, qcutoff2;
  int hp_cutoff;
} AddPathsJob;

extern boolean print_traversed_inserts;

static inline void add_path_job_alloc(AddPathsJob *job) {
  if(!seq_read_alloc(&job->r1) || !seq_read_alloc(&job->r2))
    die("Out of memory");
}

static inline void add_path_job_dealloc(AddPathsJob *job) {
  seq_read_dealloc(&job->r1);
  seq_read_dealloc(&job->r2);
}

void add_read_paths_init();
void add_read_paths_cleanup();

// Must have called add_paths_init before calling add_read_paths
void add_read_paths(const AddPathsJob *job, dBNodeBuffer *nodebuf,
                    GraphWalker *wlk, RepeatWalker *rptwlk,
                    uint64_t *insert_sizes, uint64_t *gap_sizes,
                    dBGraph *db_graph);

// Save gap size distribution
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void dump_gap_sizes(const char *base_fmt, const uint64_t *arr, size_t arrlen,
                    size_t kmer_size, boolean insert_sizes, size_t nreads);

#endif /* ADD_READ_PATH_H_ */
