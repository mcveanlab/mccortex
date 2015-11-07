#ifndef DECOMP_BUBBLE_H_
#define DECOMP_BUBBLE_H_

#include "aligned_call.h"
#include "call_file_reader.h"
#include "seq_reader.h" // genome hash

typedef struct {
  uint64_t nflank5p_unmapped, nflank5p_lowqual;
  uint64_t nflank3p_multihits, nflank3p_exact_match;
  uint64_t nflank3p_align, nflank3p_not_found, nflank3p_approx_found;
  uint64_t nflanks_overlap_too_much, nentries_well_mapped;
} DecompBubbleStats;

typedef struct DecompBubbleStruct DecompBubble;

DecompBubble* decomp_bubble_init();
void decomp_bubble_destroy(DecompBubble *bd);

void decomp_bubble_cpy_stats(DecompBubbleStats *stats, const DecompBubble *bd);

// Convert a call into an aligned call
// return 0 on success, otherwise non-zero on failure
int decomp_bubble_call(DecompBubble *bd, khash_t(ChromHash) *genome,
                       size_t kmer_size, size_t min_mapq,
                       const CallFileEntry *centry,
                       const bam1_t *mflank, const bam_hdr_t *bhdr,
                       AlignedCall *ac);

#endif /* DECOMP_BUBBLE_H_ */
