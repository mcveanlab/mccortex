#include "global.h"
#include "decomp_breakpoint.h"
#include "chrom_pos_list.h" // Parse chromosome position lists
#include "dna.h"
#include "util.h"

struct DecompBreakpointStruct
{
  ChromPosBuffer chrposbuf;
  size_t *cols, s_cols;
  DecompBreakpointStats stats;
};

DecompBreakpoint* decomp_brkpt_init()
{
  DecompBreakpoint *db = ctx_calloc(1, sizeof(*db));
  chrompos_buf_alloc(&db->chrposbuf, 32);
  return db;
}
void decomp_brkpt_destroy(DecompBreakpoint *db)
{
  chrompos_buf_dealloc(&db->chrposbuf);
  free(db->cols);
  ctx_free(db);
}

void decomp_brkpt_cpy_stats(DecompBreakpointStats *stats,
                            const DecompBreakpoint *db)
{
  memcpy(stats, &db->stats, sizeof(*stats));
}

//
// Decompose
//


/**
 * Fetch the largest match from a breakpoint call
 * @param line       input to be parsed '>seqname ... chr=...'
 * @param buf        temporary buffer
 * @param use_first  if true, return lowest offset, otherwise highest offset
 * @param flank      used to return result of largest match
 * @return true on success, false if not mapped. Calls die() on error
 */
static bool brkpnt_fetch_largest_match(const char *line, ChromPosBuffer *buf,
                                       bool use_first, ChromPosOffset *flank)
{
  char *list = strstr(line, " chr=");
  if(list == NULL) die("Cannot find flank position: %s", line);
  // Parse chr=seq0b:1-20:+:1,seq0a:2-20:+:2
  if(chrom_pos_list_parse(list+5, buf) < 0) die("Invalid positions: %s", line);
  size_t n = chrom_pos_list_get_largest(buf, use_first, flank);
  return (n > 0); // can check if n==1 if we want unique
}

// Convert a call into an aligned call
// return 0 on success, otherwise non-zero on failure
int decomp_brkpt_call(DecompBreakpoint *db,
                      ChromHash *genome, size_t nsamples,
                      const CallFileEntry *centry,
                      AlignedCall *ac)
{
  ChromPosOffset match5p, match3p;
  size_t nlines, nalleles;
  nlines = call_file_num_lines(centry);
  nalleles = (call_file_num_lines(centry) - 4) / 2;
  if(nlines < 6) die("Fewer than 6 lines: %zu", nlines);

  const char *line0 = call_file_get_line(centry, 0);
  const char *line2 = call_file_get_line(centry, 2);
  int64_t callid = call_file_get_call_id(line0, ">brkpnt.call");
  if(callid < 0) die("Bad call line: %s [%"PRIi64"]", line0, callid);

  db->stats.ncalls++;

  // Set unmapped
  ac->chrom = NULL;

  acall_resize(ac, nalleles, nsamples);

  bool good = (brkpnt_fetch_largest_match(line0, &db->chrposbuf, false, &match5p) &&
               brkpnt_fetch_largest_match(line2, &db->chrposbuf, true, &match3p));

  // Didn't map uniquely, with mismatching chromosomes or strands
  if(!good) {
    db->stats.nflanks_not_uniquely_mapped++;
    return -1;
  }
  else if(strcmp(match5p.chrom, match3p.chrom) != 0) {
    db->stats.nflanks_diff_chroms++;
    return -2;
  }
  else if(match5p.fw_strand != match3p.fw_strand) {
    db->stats.nflanks_diff_strands++;
    return -3;
  }

  // Success
  const char *flank5p = call_file_get_line(centry,1);
  const char *flank3p = call_file_get_line(centry,3);
  size_t flank5plen = call_file_line_len(centry,1);
  size_t flank3plen = call_file_line_len(centry,3);
  bool fw_strand = match5p.fw_strand;

  // sanity check
  ctx_assert(match5p.offset+chrom_pos_len(&match5p) <= flank5plen);

  // Copy required bases so flank5p, flank3p go right up to the breakpoints
  size_t flank5ptrim = flank5plen - (match5p.offset+chrom_pos_len(&match5p));
  size_t flank3ptrim = match3p.offset;

  const read_t *chrom = seq_fetch_chrom(genome, match5p.chrom);

  // ChromPosOffset.start/end are inclusive, but AlignedCall.start/end are not
  ac->start = (fw_strand ? match5p.end : match3p.end);
  ac->end   = (fw_strand ? match3p.start : match5p.start);

  if(ac->start > ac->end)
  {
    // Trim back either flank
    size_t diff, trim5p, trim3p;
    diff = ac->start - ac->end;
    trim5p = MIN2(diff, flank5plen); flank5ptrim += trim5p; diff -= trim5p;
    trim3p = MIN2(diff, flank3plen); flank3ptrim += trim3p; diff -= trim3p;
    if(diff > 0) { db->stats.nflanks_overlap_too_much++; return -4; }
    if(fw_strand) { ac->start -= trim5p; ac->end += trim3p; }
    else          { ac->start -= trim3p; ac->end += trim5p; }
  }

  // construct aligned call
  StrBuf *branch;
  const char *hdrline, *allele;
  size_t i, alen;
  int j, ncols;

  gca_resize(db->cols, db->s_cols, nsamples);
  memset(ac->gts, 0, nsamples * nalleles * sizeof(ac->gts[0]));

  for(i = 0; i < nalleles; i++)
  {
    branch = &ac->lines[i];
    hdrline = call_file_get_line(centry, 4+i*2);
    allele  = call_file_get_line(centry, 4+i*2+1);
    alen    = call_file_line_len(centry, 4+i*2+1);
    strbuf_reset(branch);
    strbuf_append_strn(branch, flank5p+flank5plen-flank5ptrim, flank5ptrim);
    strbuf_append_strn(branch, allele, alen);
    strbuf_append_strn(branch, flank3p, flank3ptrim);
    if(!fw_strand) dna_revcomp_str(branch->b, branch->b, branch->end);

    // Genotypes
    const char *colstr = strstr(hdrline,"cols=");
    if(!colstr) die("Cannot find cols=...: %s", hdrline);
    ncols = parse_list_sizes(db->cols, nsamples, colstr+strlen("cols="));
    if(ncols <= 0) die("Invalid line: %s", hdrline);
    for(j = 0; j < ncols; j++)
      ac->gts[i*nsamples + db->cols[j]] = 1;
  }

  // INFO
  StrBuf *info = &ac->info;
  strbuf_reset(info);
  strbuf_append_str(info, "BRKPNT=");
  strbuf_append_ulong(info, callid);

  // Only set chrom when we know call is aligned
  ac->chrom = chrom;
  db->stats.ncalls_mapped++;

  return 0;
}
