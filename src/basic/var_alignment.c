#include "global.h"
#include "dna.h"

// Returns msa_len if no more variants
size_t var_aln_msa_start(const char **alleles, size_t num_alleles,
                         size_t msa_len, size_t pos)
{
  size_t i;
  for(; pos < msa_len; pos++)
  {
    char c = alleles[0][pos];
    if(c == '-') return pos;
    for(i = 1; i < num_alleles && alleles[i][pos] == c; i++);
    if(i < num_alleles) return pos;
  }
  return msa_len;
}

size_t var_aln_msa_end(const char **alleles, size_t num_alleles,
                       size_t msa_len, size_t pos)
{
  size_t i; char c;
  for(pos++; pos < msa_len; pos++)
  {
    if((c = alleles[0][pos]) != '-')
    {
      for(i = 1; i < num_alleles && alleles[i][pos] == c; i++);
      if(i == num_alleles) return pos;
    }
  }
  return msa_len;
}

// Copy alignment from aln to allele, removing '-'
void var_aln_strip_allele(StrBuf *allele, const char *aln, size_t len)
{
  strbuf_reset(allele);
  size_t i;
  for(i = 0; i < len; i++) {
    if(aln[i] != '-')
      strbuf_append_char(allele, aln[i]);
  }
}

void var_aln_left_right_shift(const char *fl5p, size_t fl5plen,
                              const char *allele, size_t allelelen,
                              const char *fl3p, size_t fl3plen,
                              size_t *left, size_t *right)
{
  // no point checking futher than min_left
  size_t i, left_shift = 0, right_shift = 0;

  // Find left shift
  // Compare end of allele to end of flank5p
  while(left_shift < fl5plen &&
        allele[allelelen-left_shift-1] == fl5p[fl5plen-left_shift-1]) {
    left_shift++;
  }

  if(left_shift == allelelen)
  {
    // Can continue along 5p flank
    for(i = 0; left_shift < fl5plen &&
               fl5p[fl5plen-i-1] == fl5p[fl5plen-left_shift-1]; i++) {
      left_shift++;
    }
  }

  // Find right shift
  // Compare start of allele to start of flank3p
  while(right_shift < fl3plen && allele[right_shift] == fl3p[right_shift]) {
    right_shift++;
  }

  if(right_shift == allelelen)
  {
    // Can continue along 3p flank
    for(i = 0; right_shift < fl3plen && fl3p[i] == fl3p[right_shift]; i++)
      right_shift++;
  }

  *left = left_shift;
  *right = right_shift;
}

#ifdef VAR_ALN_DEV

//
// Decompose into variants
//

// Returns 0 if failed to find 3p flank, 1 otherwise
int var_aln_find_end_flank(const bam1_t *bam,
                           const char *flank, size_t flanklen,
                           bool is_right_fl,
                           size_t kmer_size, size_t longest_allele_bp)
{
  // check orientation, offset
  // if rev-orient, revcmp alleles and flanks
  int cigar2rlen = bam_cigar2rlen(bam->core.n_cigar, bam_get_cigar(bam));
  size_t i, endfl_missing;

  if(bam_is_rev(bam)) {
    endfl_missing = kmer_size - (lf->len-FPREFIX);
    strbuf_append_strn(&endflank, lf->buff+FPREFIX, lf->len-FPREFIX);
    strbuf_append_strn(&endflank, rf->buff+FPREFIX, endfl_missing);
  } else {
    endfl_missing = kmer_size - (rf->len-FPREFIX);
    strbuf_append_strn(&endflank, lf->buff+lf->len-endfl_missing, endfl_missing);
    strbuf_append_strn(&endflank, rf->buff+FPREFIX, rf->len-FPREFIX);
  }

  // Choose a region of the ref to search for the end flank
  // end is index after last char
  long search_start, search_end;
 
  if(bam_is_rev(bam)) {
    search_start = (long)bam->core.pos - (long)(longest_allele_bp + kmer_size + 10);
    search_end = (long)bam->core.pos + (long)kmer_size;
  } else {
    search_start = (long)(bam->core.pos + cigar2rlen) - (long)kmer_size;
    search_end = search_start + (long)(kmer_size + longest_allele_bp + kmer_size + 10);
  }

  if(search_start < 0) search_start = 0;
  if(search_end > (signed)chr->seq.end) search_end = (long)chr->seq.end;

  char *search_region = chr->seq.b+search_start;
  size_t search_len = (size_t)(search_end - search_start);

  size_t search_trim_left = 0, search_trim_right = 0;
  size_t flank_trim_left = 0, flank_trim_right = 0;

  // Attempt to find perfect match for kmer
  // temporarily null terminate ref
  char save_ref_base = search_region[search_len];
  search_region[search_len] = '\0';

  char *kmer_match = strstr(search_region, endflank.buff), *search = kmer_match;
  if(!bam_is_rev(bam) && kmer_match != NULL) {
    while((search = strstr(search+1, endflank.buff)) != NULL) {
      kmer_match = search;
    }
  }

  search_region[search_len] = save_ref_base;

  if(kmer_match != NULL)
  {
    // Found exact match
    search_trim_left = (size_t)(kmer_match - search_region);
    search_trim_right = (size_t)(&search_region[search_len]-&kmer_match[kmer_size]);

    #ifdef CTXVERBOSE
      printf("sr:%.*s\nsr:", (int)search_len, search_region);
      size_t pos = kmer_match - search_region;
      for(i = 0; i < pos; i++) fputc('.', stdout);
      fputs(endflank.buff, stdout);
      for(i = pos+kmer_size; i < (unsigned)search_len; i++) fputc('.', stdout);
      fputc('\n', stdout);
    #endif
  }
  else
  {
    // Look for approximate match
    needleman_wunsch_align2(search_region, endflank.buff, search_len, kmer_size,
                            nw_scoring_flank, nw_aligner, alignment);

    char *r1 = alignment->result_a, *r2 = alignment->result_b;

    #ifdef CTXVERBOSE
      printf("nw:%s\nnw:%s\n", r1, r2);
    #endif

    // --aa--cc-cge
    // aa--ccd-dcge

    // Find positions of first and last match
    size_t l, r, matches = 0;
    for(l = 0; l < alignment->length && r1[l] != r2[l]; l++) {
      search_trim_left += (r1[l] != '-');
      flank_trim_left += (r2[l] != '-');
    }

    for(r = alignment->length-1; r != SIZE_MAX && r1[r] != r2[r]; r--) {
      search_trim_right += (r1[r] != '-');
      flank_trim_right += (r2[r] != '-');
    }

    if(r != SIZE_MAX) {
      for(i = l; i <= r; i++) matches += (r1[i] == r2[i]);
    }

    if(matches < kmer_size / 2) {
      // flank doesn't map well
      return 0;
    }
  }

  size_t refpos, reflen;

  if(bam_is_rev(bam)) {
    reflen = search_trim_right < kmer_size ? 0 : search_trim_right - kmer_size;
    refpos = (size_t)(search_region + search_len - (reflen+kmer_size) - chr->seq.b);
    // Add flank_trim_right to the beginning of each allele
    const char *end = endflank.buff+endflank.len-flank_trim_right;
    for(i = 0; i < invcf->num_alts; i++)
      strbuf_insert(&invcf->alts[i], 0, end, flank_trim_right);
  }
  else {
    refpos = (size_t)(bam->core.pos + cigar2rlen);
    reflen = search_trim_left < kmer_size ? 0 : search_trim_left - kmer_size;
    // Append flank_trim_left onto the end of each allele
    for(i = 0; i < invcf->num_alts; i++)
      strbuf_append_strn(&invcf->alts[i], endflank.buff, flank_trim_left);
  }

  const char *ref_allele_str = chr->seq.b + refpos;
  return 1;
}

static void parse_alignment(char **alleles, size_t num_alleles, size_t msa_len,
                            const read_t *chr, size_t refpos,
                            vcf_entry_t *outvcf, FILE *fout)
{
  // Find where alleles differ
  size_t start, end = 0, j;
  size_t i, refallelelen;

  if(outvcf->lf != NULL) strbuf_reset(outvcf->lf);
  if(outvcf->rf != NULL) strbuf_reset(outvcf->rf);

  while((start = get_var_start(alleles, num_alleles, msa_len, end)) < msa_len)
  {
    // Update allele offsets
    for(j = end; j < start; j++)
      if(alleles[0][j] != '-') refpos++;

    end = get_var_end(alleles, num_alleles, msa_len, start);

    #ifdef CTXVERBOSE
      printf("start-end: %zu-%zu\n", start, end);
    #endif

    strip_allele(&outvcf->cols[VCFREF], alleles[0]+start, end-start);
    bool is_snp = (outvcf->cols[VCFREF].len == 1);
    refallelelen = outvcf->cols[VCFREF].len;

    outvcf->num_alts = num_alleles - 1;

    for(i = 0; i < outvcf->num_alts; i++) {
      strip_allele(&outvcf->alts[i], alleles[i+1]+start, end-start);
      ctx_assert(outvcf->alts[i].len == strlen(outvcf->alts[i].buff));
      is_snp &= (outvcf->alts[i].len == 1);
    }

    // vcf_info_tag_del(outvcf, "SNP");

    // REF position (still 0-based)
    long pos = (long)refpos;

    if(!is_snp)
    {
      char padding = pos > 0 ? chr->seq.b[pos-1] : 'N';
      strbuf_insert(&outvcf->cols[VCFREF], 0, &padding, 1);
      for(i = 0; i < outvcf->num_alts; i++) {
        strbuf_insert(&outvcf->alts[i], 0, &padding, 1);
      }
      pos--;
    }
    // else vcf_info_tag_add(outvcf, "SNP");

    // Collapse down matching alleles
    remove_duplicate_alleles(outvcf);

    // Set ref position
    strbuf_reset(&outvcf->cols[VCFPOS]);
    strbuf_sprintf(&outvcf->cols[VCFPOS], "%li", pos+1);

    print_entry(outvcf, fout);

    // Update ref offset
    refpos += refallelelen;
  }

  #ifdef CTXVERBOSE
    printf("//\n");
  #endif
}

// Pairwise alignment to ref
void var_aln_print_vcf(const char *seq, size_t seqlen,
                       const char *chr_name, const char *chr_seq,
                       size_t chr_start, size_t chr_end)
{
  size_t reflen = chr_end - chr_start + 1;

  if(reflen != seqlen || strcmp(chr_seq+chr_start, seq) != 0)
  {
    needleman_wunsch_align2(chr_seq+chr_start, seq, reflen, seqlen,
                            nw_scoring_allele, nw_aligner, alignment);

    char *alignments[2] = {alignment->result_a, alignment->result_b};

    var_aln_decompose(alignments, 2, alignment->length,
                      chr_name, chr_seq, chr_start, chr_end,
                      fout);
  }
}

#endif
