#ifndef DECOMP_BUBBLE_H_
#define DECOMP_BUBBLE_H_

#include "aligned_call.h"
#include "call_file_reader.h"
#include "seq_reader.h" // genome hash

typedef struct {
  uint64_t nflank5p_unmapped, nflank5p_lowqual;
  uint64_t nflank3p_multihits, nflank3p_not_found;
  uint64_t nflank3p_exact_found, nflank3p_approx_found;
  uint64_t nflanks_overlap_too_much;
  uint64_t ncalls, ncalls_mapped;
} DecompBubbleStats;

typedef struct DecompBubbleStruct DecompBubble;

DecompBubble* decomp_bubble_init();
void decomp_bubble_destroy(DecompBubble *db);

void decomp_bubble_cpy_stats(DecompBubbleStats *stats, const DecompBubble *db);
scoring_t* decomp_bubble_get_scoring(DecompBubble *db);

// Convert a call into an aligned call
// return 0 on success, otherwise non-zero on failure
int decomp_bubble_call(DecompBubble *db, ChromHash *genome,
                       size_t kmer_size, size_t min_mapq,
                       const CallFileEntry *centry,
                       const bam1_t *mflank, const bam_hdr_t *bhdr,
                       AlignedCall *ac);

#endif /* DECOMP_BUBBLE_H_ */
