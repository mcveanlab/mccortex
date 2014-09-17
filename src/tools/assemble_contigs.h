#ifndef ASSEMBLE_CONTIGS_H_
#define ASSEMBLE_CONTIGS_H_

#include "db_graph.h"
#include "graph_walker.h"
#include "common_buffers.h"

#include "seq_file.h"

#define AC_MAX_PATHS 5

typedef struct
{
  uint64_t num_contigs, total_len, total_junc;
  uint64_t contigs_outdegree[5];
  uint64_t paths_held[AC_MAX_PATHS]; // Paths held when we stopped
  uint64_t paths_new[AC_MAX_PATHS];
  uint64_t paths_cntr[AC_MAX_PATHS];
  uint64_t paths_held_max, paths_new_max, paths_cntr_max;
  uint64_t grphwlk_steps[GRPHWLK_NUM_STATES]; // states in graph_walker.h
  // uint64_t *lengths, *junctns, capacity; // length, njuncs for each contig
  SizeBuffer lengths, junctns;
  double max_junc_density;
  uint64_t num_reseed_abort; // aborted - already visited seed
  uint64_t num_seeds_not_found; // seed contig didn't have any matching kmers
} AssembleContigStats;

void assemble_contigs_stats_init(AssembleContigStats *stats);
void assemble_contigs_stats_destroy(AssembleContigStats *stats);
void assemble_contigs_stats_print(const AssembleContigStats *stats);
void assemble_contigs_stats_merge(AssembleContigStats *dst,
                                  const AssembleContigStats *src);

/*!
  @param contig_limit Stop after printing this many contigs, if zero no limit
  @param visited Bit array to store visited nodes in. If not NULL, do not use a
                 seed that has previously been visited. We do not clear this
                 array.
 */
void assemble_contigs(size_t nthreads,
                      seq_file_t **seed_files, size_t num_seed_files,
                      size_t contig_limit, uint8_t *visited,
                      FILE *fout, const char *out_path,
                      AssembleContigStats *stats,
                      const dBGraph *db_graph, size_t colour);

#endif /* ASSEMBLE_CONTIGS_H_ */
