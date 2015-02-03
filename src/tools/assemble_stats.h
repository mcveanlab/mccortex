#ifndef ASSEMBLE_STATS_H_
#define ASSEMBLE_STATS_H_

#include "db_graph.h"
#include "common_buffers.h"
#include "graph_walker.h"

enum AssemStopCause {
  ASSEM_STOP_UNKNOWN        = 0,
  ASSEM_STOP_NOCOVG         = 1, /* Fail: no choices */
  ASSEM_STOP_NOCOLCOVG      = 2, /* Fail: fork in pop but no choices in colour */
  ASSEM_STOP_NOPATHS        = 3, /* Fail: fork in colour, no paths */
  ASSEM_STOP_SPLIT_PATHS    = 4, /* Fail: fork in colour, paths split */
  ASSEM_STOP_MISSING_PATHS  = 5, /* Fail: fork in colour, missing info */
  ASSEM_STOP_CYCLE          = 6,
  ASSEM_STOP_LOW_STEP_CONF  = 7,
  ASSEM_STOP_LOW_CUMUL_CONF = 8
};

#define ASSEM_NUM_STOPS 9

#define ASSEM_STOP_UNKNOWN_STR        "StopUnknown"
#define ASSEM_STOP_NOCOVG_STR         "StopNoCovg"
#define ASSEM_STOP_NOCOLCOVG_STR      "StopPopForkNoColCovg"
#define ASSEM_STOP_NOPATHS_STR        "StopForkNoPaths"
#define ASSEM_STOP_SPLIT_PATHS_STR    "StopForkInPaths"
#define ASSEM_STOP_MISSING_PATHS_STR  "StopMissingPaths"
#define ASSEM_STOP_CYCLE_STR          "StopHitLoop"
#define ASSEM_STOP_LOW_STEP_CONF_STR  "StopLowStepConfidence"
#define ASSEM_STOP_LOW_CUMUL_CONF_STR "StopLowCumulativeConfidence"

enum AssemStopCause graphstep2assem(enum GraphStepStatus step, bool hit_cycle,
                                    bool low_step_confid, bool low_cumul_confid);

char* assem2str(enum AssemStopCause assem, char *str, size_t size);

#define AC_MAX_PATHS 5

typedef struct
{
  uint64_t num_contigs, total_len, total_junc;
  uint64_t num_rnd_seeds, num_gpath_seeds;
  uint64_t contigs_outdegree[5];
  uint64_t paths_held[AC_MAX_PATHS]; // Paths held when we stopped
  uint64_t paths_new[AC_MAX_PATHS];
  uint64_t paths_cntr[AC_MAX_PATHS];
  uint64_t paths_held_max, paths_new_max, paths_cntr_max;
  uint64_t grphwlk_steps[GRPHWLK_NUM_STATES]; // states in graph_walker.h
  uint64_t stop_causes[ASSEM_NUM_STOPS]; // ASSEM_STOP_* defined above
  SizeBuffer lengths, junctns;
  double max_junc_density;
  uint64_t num_contigs_from_seed_kmers;
  uint64_t num_contigs_from_seed_paths;
  uint64_t num_reseed_abort; // aborted - already visited seed
  uint64_t num_seeds_not_found; // seed contig didn't have any matching kmers
} AssembleContigStats;

// Results from a single contig
struct ContigStats {
  size_t num_nodes, num_junc;
  size_t wlk_steps[GRPHWLK_NUM_STATES];
  size_t paths_held[2], paths_cntr[2];
  uint8_t stop_causes[2]; // one of ASSEM_STOP_* for each direction
  size_t max_step_gap[2];
  double gap_conf[2];
  uint8_t outdegree_rv, outdegree_fw;
  bool seed_kmer, seed_path;
  size_t num_seed_kmers; // number of kmers used to seed
};

void assemble_contigs_stats_init(AssembleContigStats *stats);
void assemble_contigs_stats_destroy(AssembleContigStats *stats);
void assemble_contigs_stats_print(const AssembleContigStats *stats);
void assemble_contigs_stats_add(AssembleContigStats *stats,
                                const struct ContigStats *s);
void assemble_contigs_stats_merge(AssembleContigStats *dst,
                                  const AssembleContigStats *src);

#endif /* ASSEMBLE_STATS_H_ */
