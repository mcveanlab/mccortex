#ifndef BUILD_GRAPH_H_
#define BUILD_GRAPH_H_

//
// Multithreaded graph building
// see generate_paths.h for a path threading comparison
//

#include "msg-pool/msgpool.h"
#include "seq_file.h"

#include "cortex_types.h"
#include "db_graph.h"
#include "seq_reader.h"
#include "async_read_io.h"
#include "loading_stats.h"

typedef struct
{
  AsyncIOReadTask files;
  const uint8_t fq_cutoff, hp_cutoff;
  const ReadMateDir matedir;
  const Colour colour;
  const bool remove_pcr_dups;

  // Stats are written to here
  LoadingStats stats;

  // used internally
  size_t idx;
} BuildGraphTask;

void build_graph_task_print(const BuildGraphTask *task);
void build_graph_task_print_stats(const BuildGraphTask *task);

// Threadsafe graph construction
// Beware: this function does not update ginfo
void build_graph_from_reads_mt(read_t *r1, read_t *r2,
                               uint8_t fq_offset1, uint8_t fq_offset2,
                               uint8_t fq_cutoff, uint8_t hp_cutoff,
                               bool remove_pcr_dups, ReadMateDir matedir,
                               LoadingStats *stats, size_t colour,
                               dBGraph *db_graph);

// One thread used per input file, num_build_threads used to add reads to graph
// Updates ginfo
void build_graph(dBGraph *db_graph, BuildGraphTask *files,
                 size_t num_files, size_t num_build_threads);

// Threadsafe
// Sequence must be entirely ACGT and len >= kmer_size
// Returns number of novel kmers loaded
size_t build_graph_from_str_mt(dBGraph *db_graph, size_t colour,
                               const char *seq, size_t len);

#endif /* BUILD_GRAPH_H_ */
