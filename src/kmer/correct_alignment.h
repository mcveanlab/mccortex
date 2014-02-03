#ifndef CORRECTED_ALIGNMENT_H_
#define CORRECTED_ALIGNMENT_H_

#include "cortex_types.h"
#include "db_graph.h"
#include "db_node.h"
#include "db_alignment.h"
#include "graph_walker.h"
#include "repeat_walker.h"

typedef struct
{
  Colour ctpcol, ctxcol;
  uint32_t ins_gap_min, ins_gap_max;
  boolean one_way_gap_traverse; // set to false for more error prone algo
  uint32_t max_context, gap_wiggle;
  float gap_variance; // permitted gap size = X*gap_variance + gap_wiggle
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
  uint64_t *gap_ins_histgrm, *gap_err_histgrm, histgrm_len;
} CorrectAlnWorker;

size_t correct_aln_worker_est_mem(const dBGraph *graph);

void correct_aln_worker_alloc(CorrectAlnWorker *wrkr, const dBGraph *db_graph);
void correct_aln_worker_dealloc(CorrectAlnWorker *wrkr);

void correct_alignment_init(CorrectAlnWorker *wrkr, const dBAlignment *aln,
                            CorrectAlnParam params);

// DEV: add option to extend
// Returns NULL if end of alignment
dBNodeBuffer* correct_alignment_nxt(CorrectAlnWorker *wrkr);


// Get alignment coords of contig
size_t correct_alignment_get_strtidx(CorrectAlnWorker *wrkr);
size_t correct_alignment_get_endidx(CorrectAlnWorker *wrkr);

uint64_t* correct_alignment_get_errhist(CorrectAlnWorker *wrkr, size_t *n);
uint64_t* correct_alignment_get_inshist(CorrectAlnWorker *wrkr, size_t *n);

// copy to dst histrograms, zero src histograms
void correct_alignment_merge_hists(CorrectAlnWorker *dst, CorrectAlnWorker *src);

#endif /* CORRECTED_ALIGNMENT_H_ */
