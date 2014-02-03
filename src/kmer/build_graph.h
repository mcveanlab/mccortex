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
  const boolean remove_dups_se, remove_dups_pe;
  const ReadMateDir matedir;

  // Stats are written to here
  LoadingStats stats;

  // used internally
  size_t idx;
} BuildGraphTask;

// One thread used per input file, num_build_threads used to add reads to graph
void build_graph(dBGraph *db_graph, BuildGraphTask *files,
                 size_t num_files, size_t num_build_threads);

// Threadsafe
// Sequence must be entirely ACGT and len >= kmer_size
void build_graph_from_str_mt(dBGraph *db_graph, size_t colour,
                             const char *seq, size_t len);

#endif /* BUILD_GRAPH_H_ */
