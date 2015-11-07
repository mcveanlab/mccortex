#ifndef DECOMP_BREAKPOINT_H_
#define DECOMP_BREAKPOINT_H_

#include "aligned_call.h"
#include "call_file_reader.h"
#include "seq_reader.h" // genome hash

typedef struct {
  uint64_t nflanks_not_uniquely_mapped, nflanks_diff_chroms;
  uint64_t nflanks_diff_strands, nflanks_overlap_too_much;
} DecompBreakpointStats;

typedef struct DecompBreakpointStruct DecompBreakpoint;

DecompBreakpoint* decomp_brkpt_init();
void decomp_brkpt_destroy(DecompBreakpoint *bd);

void decomp_brkpt_cpy_stats(DecompBreakpointStats *stats,
                            const DecompBreakpoint *bd);

// Convert a call into an aligned call
// return 0 on success, otherwise non-zero on failure
int decomp_brkpt_call(DecompBreakpoint *db,
                      khash_t(ChromHash) *genome, size_t nsamples,
                      const CallFileEntry *centry,
                      AlignedCall *ac);

#endif /* DECOMP_BREAKPOINT_H_ */
