#ifndef CORRECTED_ALIGNMENT_H_
#define CORRECTED_ALIGNMENT_H_

#include "cortex_types.h"
#include "db_graph.h"
#include "db_node.h"
#include "db_alignment.h"
#include "graph_walker.h"
#include "repeat_walker.h"
#include "correct_aln_stats.h"

#define DEFAULT_CRTALN_FRAGLEN_MIN 0
#define DEFAULT_CRTALN_FRAGLEN_MAX 1000

// seq gap of N bases can be filled by MAX2(0, N±(N*GAP_VARIANCE+GAP_WIGGLE))
#define DEFAULT_CRTALN_GAP_VARIANCE 0.1
#define DEFAULT_CRTALN_GAP_WIGGLE 5

#define DEFAULT_CRTALN_MAX_CONTEXT 200

#define CORRECT_PARAMS_DEFAULT {.ctpcol = 0, .ctxcol = 0,                      \
                                .frag_len_min = DEFAULT_CRTALN_FRAGLEN_MIN,    \
                                .frag_len_max = DEFAULT_CRTALN_FRAGLEN_MAX,    \
                                .one_way_gap_traverse = true,                  \
                                .use_end_check = true,                         \
                                .max_context = DEFAULT_CRTALN_MAX_CONTEXT,     \
                                .gap_variance = DEFAULT_CRTALN_GAP_VARIANCE,   \
                                .gap_wiggle = DEFAULT_CRTALN_GAP_WIGGLE}

typedef struct
{
  Colour ctpcol, ctxcol;
  // uint32_t ins_gap_min, ins_gap_max;
  uint32_t frag_len_min, frag_len_max; // For PE reads
  uint32_t max_context; // how many kmers to use either side of a gap
  float gap_wiggle, gap_variance; // permitted gap size = X*gap_variance + gap_wiggle
  bool one_way_gap_traverse; // set to false for more error prone algo
  bool use_end_check; // check paths match remaining contig
} CorrectAlnParam;

typedef struct
{
  const dBGraph *const db_graph;
  GraphWalker wlk, wlk2;
  RepeatWalker rptwlk, rptwlk2;

  // Alignment
  const dBAlignment *aln;
  CorrectAlnParam params;

  // Current State
  //
  // start_idx  gap_idx  end_idx
  // v..........v........v
  //
  // with no gap:
  // start_idx           gap_idx|end_idx
  // v...................v
  // return alignment from [start_idx..gap_idx]
  size_t start_idx, gap_idx, end_idx;
  size_t prev_start_idx;

  // contig with gaps filled
  // we use revcontig when walking backwards
  dBNodeBuffer contig, revcontig;

  // Statistics on gap traversal
  LoadingStats stats;
  CorrectAlnStats gapstats;
} CorrectAlnWorker;

// Global variables to specify if we should print output - used for debugging only
// These are used in generate_paths.c
extern bool gen_paths_print_contigs, gen_paths_print_paths, gen_paths_print_reads;


size_t correct_aln_worker_est_mem(const dBGraph *graph);

void correct_aln_worker_alloc(CorrectAlnWorker *wrkr, const dBGraph *db_graph);
void correct_aln_worker_dealloc(CorrectAlnWorker *wrkr);

// Merge stats into dst and reset src
void correct_aln_merge_stats(CorrectAlnWorker *restrict dst,
                             CorrectAlnWorker *restrict src);

void correct_alignment_init(CorrectAlnWorker *wrkr, const dBAlignment *aln,
                            CorrectAlnParam params);

// DEV: add option to extend by n bases for read correction?
// Returns NULL if end of alignment
dBNodeBuffer* correct_alignment_nxt(CorrectAlnWorker *wrkr);

// Get alignment coords of contig
// Called after correct_alignment_nxt()
size_t correct_alignment_get_strtidx(CorrectAlnWorker *wrkr);
size_t correct_alignment_get_endidx(CorrectAlnWorker *wrkr);

#endif /* CORRECTED_ALIGNMENT_H_ */
