#ifndef ASSEMBLE_STATS_H_
#define ASSEMBLE_STATS_H_

#include "db_graph.h"
#include "common_buffers.h"
#include "graph_walker.h"

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
  uint64_t num_cycles; // instances where the repeat_walker stopped us
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

#endif /* ASSEMBLE_STATS_H_ */
