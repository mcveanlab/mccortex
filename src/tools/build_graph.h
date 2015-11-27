#ifndef BUILD_GRAPH_H_
#define BUILD_GRAPH_H_

//
// Multithreaded graph building
// see generate_paths.h for a path threading comparison
//

#include "msg-pool/msgpool.h"
#include "seq_file/seq_file.h"

#include "cortex_types.h"
#include "db_graph.h"
#include "seq_reader.h"
#include "async_read_io.h"
#include "seq_loading_stats.h"

typedef struct
{
  uint8_t fq_cutoff, hp_cutoff;
  ReadMateDir matedir;
  Colour colour;
  bool remove_pcr_dups, must_exist_in_graph;
} SeqLoadingPrefs;

typedef struct
{
  AsyncIOInput files;
  SeqLoadingPrefs prefs;

  // Stats are written to here
  SeqLoadingStats stats;

  // used internally
  size_t idx;
} BuildGraphTask;

#define SEQ_LOADING_PREFS_INIT (SeqLoadingPrefs){.fq_cutoff = 0, \
                                                 .hp_cutoff = 0, \
                                                 .matedir = READPAIR_FR, \
                                                 .colour = 0, \
                                                 .remove_pcr_dups = false}

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(build_graph_task_buf, BuildGraphTaskBuffer, BuildGraphTask);

static inline void build_graph_task_destroy(BuildGraphTask *task)
{
  asyncio_task_close(&task->files);
}

void build_graph_task_print(const BuildGraphTask *task);
void build_graph_task_print_stats(const BuildGraphTask *task);

// Threadsafe graph construction
// Beware: this function does not update ginfo
void build_graph_from_reads_mt(read_t *r1, read_t *r2,
                               uint8_t fq_offset1, uint8_t fq_offset2,
                               const SeqLoadingPrefs *prefs,
                               SeqLoadingStats *stats,
                               dBGraph *db_graph);

// One thread used per input file, num_build_threads used to add reads to graph
// Updates ginfo
void build_graph(dBGraph *db_graph, BuildGraphTask *files,
                 size_t num_files, size_t num_build_threads);

// One thread used per input file, num_build_threads used to add reads to graph
// Updates ginfo
void build_graph_from_seq(dBGraph *db_graph, seq_file_t **files,
                          size_t num_files, size_t num_build_threads,
                          size_t colour);

// Threadsafe
// Sequence must be entirely ACGT and len >= kmer_size
// Returns number of non-novel kmers seen
size_t build_graph_from_str_mt(dBGraph *db_graph, size_t colour,
                               const char *seq, size_t len,
                               bool must_exist_in_graph);

#endif /* BUILD_GRAPH_H_ */
