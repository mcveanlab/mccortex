#ifndef BUILD_GRAPH_H_
#define BUILD_GRAPH_H_

//
// Multithreaded graph building
// see generate_paths.h for a path threading comparison
//

#include "msgpool.h"
#include "seq_file.h"

#include "cortex_types.h"
#include "db_graph.h"
#include "seq_reader.h"
#include "loading_stats.h"

typedef struct
{
  seq_file_t *const file1, *const file2;
  const Colour colour;
  // FASTQ qual offset, cutoff and homopolymer cutoff
  const uint8_t fq_offset, fq_cutoff, hp_cutoff;
  const bool remove_pcr_dups;
  const ReadMateDir matedir;

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
