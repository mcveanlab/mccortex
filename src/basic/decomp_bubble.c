#include "global.h"
#include "decomp_bubble.h"
#include "util.h"
#include "dna.h"

struct DecompBubbleStruct
{
  nw_aligner_t *nw_aligner;
  scoring_t *scoring;
  alignment_t *aln;
  StrBuf flank3pbuf;
  // stats
  DecompBubbleStats stats;
};

DecompBubble* decomp_bubble_init()
{
  DecompBubble *db = ctx_calloc(1, sizeof(DecompBubble));
  db->nw_aligner = needleman_wunsch_new();
  db->aln = alignment_create(1024);
  db->scoring = ctx_calloc(1, sizeof(db->scoring[0]));
  scoring_system_default(db->scoring);
  strbuf_alloc(&db->flank3pbuf, 1024);
  return db;
}

void decomp_bubble_destroy(DecompBubble *db)
{
  alignment_free(db->aln);
  needleman_wunsch_free(db->nw_aligner);
  ctx_free(db->scoring);
  strbuf_dealloc(&db->flank3pbuf);
  ctx_free(db);
}

void decomp_bubble_cpy_stats(DecompBubbleStats *stats, const DecompBubble *db)
{
  memcpy(stats, &db->stats, sizeof(*stats));
}

scoring_t* decomp_bubble_get_scoring(DecompBubble *db)
{
  return db->scoring;
}

//
// Decompose into aligned call
//

// returns number of bases borrowed from 5' flank
static uint32_t bubble_get_end_kmer(const char *flank5p, size_t flank5p_len,
                                    const char *flank3p, size_t flank3p_len,
                                    char *endkmer, size_t ksize)
{
  // 3p flank may not be long enough to give kmer bases
  size_t flank3pcpy = MIN2(ksize, flank3p_len);
  size_t flank5pcpy = ksize - flank3pcpy; // Make up remaining sequence
  ctx_assert(flank5pcpy <= flank5p_len);

  memcpy(endkmer,            flank5p, flank5pcpy);
  memcpy(endkmer+flank5pcpy, flank3p, flank3pcpy);
  endkmer[ksize] = '\0';

  return flank5pcpy;
}


// returns true if 3' flank aligns well, otherwise false
// *pos is set to the ref position of the last base BEFORE the flank
static bool align_flank3p(const bam1_t *mflank,
                          uint32_t cigar2rlen,
                          uint32_t max_alen, uint32_t kmer_size,
                          const char *flank5p, uint32_t flank5plen,
                          const char *flank3p, uint32_t flank3plen,
                          const char *chr, uint32_t chrlen,
                          int32_t *pos, uint32_t *flank3ptrim,
                          DecompBubble *db)
{
  // Find 3p flank position using search for first kmer
  // Determine search space
  // Choose a region of the ref to search for the end flank
  // end is index after last char
  int32_t search_start, search_end;
  bool fw_strand = !bam_is_rev(mflank);

  char endkmer[200];
  uint32_t flank5pcpy = bubble_get_end_kmer(flank5p, flank5plen,
                                            flank3p, flank3plen,
                                            endkmer, kmer_size);
  if(!fw_strand) dna_revcomp_str(endkmer, endkmer, kmer_size);

  if(fw_strand) {
    search_start = mflank->core.pos;
    search_end = mflank->core.pos + cigar2rlen + max_alen + kmer_size*2 + 10;
  } else {
    search_start = mflank->core.pos - (max_alen + kmer_size*2 + 10);
    search_end = mflank->core.pos + cigar2rlen;
  }

  search_start = MAX2(search_start, 0);
  search_end   = MIN2(search_end, (int32_t)chrlen);

  const char *search_region = chr + search_start;
  size_t search_len = (size_t)(search_end - search_start);

  // Now do search with kmer
  // Attempt to find perfect match for kmer within search region

  // Search, if there is more than one match -> abandon
  const char *kmer_match = ctx_strnstr(search_region, endkmer, search_len);

  if(kmer_match != NULL)
  {
    // Check for multiple hits
    size_t kstart = kmer_match - chr;
    if(ctx_strnstr(kmer_match+1, endkmer, search_end - kstart - 1) != NULL) {
      db->stats.nflank3p_multihits++;
      return false;
    }
    *flank3ptrim = 0;
    *pos = fw_strand ? kstart - 1 + flank5pcpy : kstart + kmer_size - flank5pcpy;
    db->stats.nflank3p_exact_found++;
    return true;
  }

  // Look for approximate match
  needleman_wunsch_align2(search_region, endkmer, search_len, kmer_size,
                          db->scoring, db->nw_aligner, db->aln);
  const char *ref = db->aln->result_a, *alt = db->aln->result_b;
  // --aa--dd-cge
  // bb--ccd-ecge

  // Find positions of first and last match
  uint32_t i, matches = 0;
  uint32_t ref_skip = 0, alt_skip = 0;

  if(fw_strand) {
    for(i = 0; i < db->aln->length; i++) {
      if(ref[i] == alt[i] && alt_skip >= flank5pcpy) break;
      ref_skip += (ref[i] != '-');
      alt_skip += (alt[i] != '-');
    }
  } else {
    for(i = db->aln->length-1; ; i--) {
      if(ref[i] == alt[i] && alt_skip >= flank5pcpy) break;
      ref_skip += (ref[i] != '-');
      alt_skip += (alt[i] != '-');
      if(i == 0) break;
    }
  }

  // Count matches
  for(i = 0; i < db->aln->length; i++) matches += (ref[i] == alt[i]);

  if(matches < kmer_size / 2 || alt_skip < flank5pcpy)
  {
    // flank doesn't map well
    db->stats.nflank3p_not_found++;
    return false;
  }

  db->stats.nflank3p_approx_found++;

  *flank3ptrim = alt_skip - flank5pcpy;
  *pos = fw_strand ? search_start + ref_skip - 1
                   : search_end - ref_skip;

  return true;
}

const uint32_t qclip = (1<<BAM_CINS)|(1<<BAM_CSOFT_CLIP)|(1<<BAM_CHARD_CLIP);

#define cigar_is_qclip(x) \
  ((qclip >> bam_cigar_op(x)) & 1)

// Return sum of bases on right of alignment with:
// * hard masked (H)
// * soft masked (S)
// * inserted bases relative to ref (I)
static inline uint32_t bam_get_end_trim(const uint32_t **_cigar, uint32_t *_n,
                                        bool fw_strand)
{
  uint32_t i, idx, bases = 0, trim = 0, n = *_n;
  const uint32_t *cigar = *_cigar;

  for(i = 0; i < n-1; i++) {
    idx = fw_strand ? n-1-i : i;
    if(cigar_is_qclip(cigar[idx])) {
      bases += bam_cigar_oplen(cigar[idx]);
      trim++;
    }
    else break;
  }

  if(!fw_strand) (*_cigar) += trim;
  (*_n) -= trim;

  return bases;
}

// Consume at least nref bases from alignment
static inline void bam_cigar_consume_ref(const uint32_t *cigar, uint32_t n,
                                         uint32_t nref, bool fw,
                                         uint32_t *qdrop, uint32_t *rdrop)
{
  ctx_assert(n > 0);
  uint32_t i, idx, len, type;
  *qdrop = *rdrop = 0;
  for(i = 0; i < n-1 && *rdrop < nref; i++)
  {
    idx = fw ? i : n-1-i;
    len = bam_cigar_oplen(cigar[idx]);
    type = bam_cigar_type(cigar[idx]); // bit 0: query, bit 1: ref
    if(type & 1) *qdrop += len;
    if(type & 2) *rdrop += len;
  }
}


// Trim up to k-1 bases from the end of bubble paths and copy to 3p flank
// return number of bases trimmed off the end of each allele
static size_t bubble_get_flank3p(const CallFileEntry *centry,
                                 size_t kmer_size, StrBuf *flank3pbuf)
{
  size_t i, trimlen, min_allele_len = SIZE_MAX;
  size_t nlines = call_file_num_lines(centry);

  for(i = 5; i < nlines; i += 2)
    min_allele_len = MIN2(call_file_line_len(centry,i), min_allele_len);

  trimlen = MIN2(min_allele_len, kmer_size-1);

  // Copy flank 3p
  const char *trimstr = call_file_get_line(centry, 5) +
                        call_file_line_len(centry, 5) - trimlen;
  strbuf_reset(flank3pbuf);
  strbuf_append_strn(flank3pbuf, trimstr, trimlen);
  strbuf_append_strn(flank3pbuf, call_file_get_line(centry, 3),
                                 call_file_line_len(centry, 3));

  return trimlen;
}

// Convert a call into an aligned call
// return 0 on success, otherwise non-zero on failure
int decomp_bubble_call(DecompBubble *db, ChromHash *genome,
                       size_t kmer_size, size_t min_mapq,
                       const CallFileEntry *centry,
                       const bam1_t *mflank, const bam_hdr_t *bhdr,
                       AlignedCall *ac)
{
  size_t nlines, nalleles;
  nlines = call_file_num_lines(centry);
  nalleles = (call_file_num_lines(centry) - 4) / 2;
  if(nlines < 6) die("Fewer than 6 lines: %zu", nlines);

  const char *line0 = call_file_get_line(centry, 0);
  int64_t callid = call_file_get_call_id(line0, ">bubble.call");
  if(callid < 0) die("Bad call line: %s [%"PRIi64"]", line0, callid);

  // Check sample name / flank mapping qname match
  const char *bname = bam_get_qname(mflank);
  const char *hdrline = call_file_get_line(centry, 0);
  if(seq_read_names_cmp(hdrline+1, bname) != 0) // hdrline[0] == '>'
    die("SAM/BAM and call entries mismatch '%s' vs '%s'", hdrline, bname);

  db->stats.ncalls++;

  // Set unmapped
  ac->chrom = NULL;

  // bubble callers don't output any sample info
  acall_resize(ac, nalleles, 0);

  // Check mapping
  if(mflank->core.flag & BAM_FUNMAP){ db->stats.nflank5p_unmapped++; return -1; }
  if(mflank->core.qual < min_mapq)  { db->stats.nflank5p_lowqual++;  return -2; }
  bool fw_strand = !bam_is_rev(mflank);
  const uint32_t *cigar = bam_get_cigar(mflank);
  uint32_t n_cigar = mflank->core.n_cigar;

  // Trim off 5p and add to alleles
  uint32_t flank5ptrim = bam_get_end_trim(&cigar, &n_cigar, fw_strand);
  const char *flank5p = call_file_get_line(centry, 1);
  uint32_t flank5plen = call_file_line_len(centry, 1);
  int32_t cigar2rlen = bam_cigar2rlen(n_cigar, cigar);
  int32_t flank5ppos = fw_strand ? mflank->core.pos + cigar2rlen
                                 : mflank->core.pos-1;

  // Find chromosome
  const char *chrom_name = bhdr->target_name[mflank->core.tid];
  const read_t *chrom = chrom_hash_fetch(genome, chrom_name);

  // Construct complete 3p flank
  // alleletrim is the number of bases taken from the right end of alleles
  uint32_t alleletrim = bubble_get_flank3p(centry, kmer_size, &db->flank3pbuf);
  // printf("alleletrim: %zu\n", alleletrim);
  uint32_t max_alen = call_file_max_allele_len(centry) - alleletrim;
  const char *flank3p = db->flank3pbuf.b;
  uint32_t flank3plen = db->flank3pbuf.end;
  // flank3plen may be as short as 1bp

  // flank5ppos, flank3ppos are the pos on ref of the base inside the flank.
  // FORWARD: <flank5p>X...Y<flank3p>   X=flank5ppos
  // REVERSE: <flank5p>Y...X<flank5p>   Y=flank3ppos
  int32_t flank3ppos = 0;
  uint32_t flank3ptrim = 0;

  // find flank and set flank3ptrim, flank3ppos
  bool mapped3p = align_flank3p(mflank, cigar2rlen, max_alen, kmer_size,
                                flank5p, flank5plen, flank3p, flank3plen,
                                chrom->seq.b, chrom->seq.end,
                                &flank3ppos, &flank3ptrim, db);

  if(!mapped3p) return -3;

  ac->start = (fw_strand ? flank5ppos : flank3ppos);
  ac->end   = (fw_strand ? flank3ppos : flank5ppos) + 1;

  uint32_t qdrop = 0, rdrop = 0, rdiff = ac->start - ac->end;

  if(ac->start > ac->end)
  {
    // Trim back 5' flank

    // if we remove diff bases from flank5p, how many ref bases do we move?
    bam_cigar_consume_ref(cigar, n_cigar, rdiff, fw_strand, &qdrop, &rdrop);

    if(rdrop < rdiff) { db->stats.nflanks_overlap_too_much++; return -4; }

    flank5ptrim += qdrop;
    flank5ppos = fw_strand ? flank5ppos - rdrop : flank5ppos + rdrop;

    ac->start = (fw_strand ? flank5ppos : flank3ppos);
    ac->end   = (fw_strand ? flank3ppos : flank5ppos) + 1;
  }

  ctx_assert(flank3ptrim <= flank3plen);
  ctx_assert(flank5ptrim <= flank5plen);

  StrBuf *branch;
  const char *allele;
  size_t i, alen;

  // Construct aligned call branches
  for(i = 0; i < nalleles; i++)
  {
    branch = &ac->lines[i];
    allele = call_file_get_line(centry, 5+i*2);
    alen = call_file_line_len(centry, 5+i*2);
    strbuf_reset(branch);
    strbuf_append_strn(branch, flank5p+flank5plen-flank5ptrim, flank5ptrim);
    strbuf_append_strn(branch, allele, alen-alleletrim);
    strbuf_append_strn(branch, flank3p, flank3ptrim);
    if(!fw_strand) dna_revcomp_str(branch->b, branch->b, branch->end);
  }

  // INFO
  StrBuf *info = &ac->info;
  strbuf_reset(info);
  strbuf_append_str(info, "BUBBLE=");
  strbuf_append_ulong(info, callid);

  // Only set chrom when we know call is aligned
  ac->chrom = chrom;
  db->stats.ncalls_mapped++;

  return 0;
}

