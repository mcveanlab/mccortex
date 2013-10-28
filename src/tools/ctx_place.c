#include "global.h"

#include <time.h>

// #include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h> // beta distribution
#include <gsl/gsl_sf_gamma.h> // gamma and beta functions

// #include "vcf.h" // not using htslib
#include "sam.h"
#include "khash.h"
#include "seq_file.h"
#include "string_buffer.h"
#include "needleman_wunsch.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "vcf_parsing.h"
#include "binary_kmer.h"
#include "file_reader.h"
#include "delta_array.h"
#include "supernode.h"

static const char usage[] =
"usage: "CMD" place [options] <calls.vcf> <calls.sam> <ref1.fa ...>\n"
"  Align calls to a reference genome.\n"
"  Options:\n"
"    --minmapq <mapq>    Flank must map with MAPQ >= <mapq> [default: 30]\n"
"    --maxbr <num>       Don't process bubbles with > <num> branches [default: 6]\n"
"    --genome <bases>    Organism genome size (e.g. 3.1G for humans)\n"
"    --ploidy <sample>:<chr>:<ploidy> [default: *:*:2]\n"
"      <sample> and <chr> can be * for all or a comma-separated list\n";

typedef struct {
  read_t r;
  uint32_t index;
} chrom_t;

// Reference genome
// Hash map of chromosome name -> sequence
KHASH_MAP_INIT_STR(ghash, chrom_t*)
khash_t(ghash) *genome;
chrom_t *chroms;
uint32_t num_chroms = 0, chrom_capacity;

// Flank mapping
bam_hdr_t *bam_header;

// VCF header info
KHASH_MAP_INIT_STR(samplehash, uint32_t)
khash_t(samplehash) *sample_indx;
uint32_t num_samples, samples_capacity;
char **sample_names;
uint64_t *sample_total_seq;

// 2d array of [chrom][sample] = ploidy[chrom*num_samples + sample]
uint32_t *ploidy;

// nw alignment
nw_aligner_t *nw_aligner;
alignment_t *alignment;
scoring_t scoring;

// Multiple sequence alignment (msa)
size_t msa_capacity = 1024, num_msa_alleles = 16;
char *msa_wrking, **msa_alleles;

// Temporary memory
StrBuf *endflank;

// Filtering parameters
uint32_t min_mapq = 30, max_branches = 6;
size_t genome_size = 0;

// VCF printing
uint32_t var_print_num = 0;

static void print_entry(vcf_entry_t *vcfentry, FILE *out)
{
  // Set var num
  strbuf_reset(&vcfentry->cols[VCFID]);
  strbuf_sprintf(&vcfentry->cols[VCFID], "var%u", var_print_num++);
  vcf_entry_print(vcfentry, out, num_samples);
}

static int get_var_start(char **alleles, int num_alleles, int msa_len, int pos)
{
  int i;
  for(; pos < msa_len; pos++)
  {
    char c = alleles[0][pos];
    if(c == '-') return pos;
    for(i = 1; i < num_alleles && alleles[i][pos] == c; i++);
    if(i < num_alleles) return pos;
  }
  return -1;
}

static int get_var_end(char **alleles, int num_alleles, int msa_len, int pos)
{
  int i;
  char c;
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

static void strip_allele(StrBuf *allele, const char* input, size_t len)
{
  strbuf_reset(allele);
  size_t i;
  for(i = 0; i < len; i++) {
    if(input[i] != '-') strbuf_append_char(allele, input[i]);
  }
}

// Do multiple sequence alignment
// Returns length, stores alleles (not null always terminated) in alleles
static uint32_t do_msa(char **seqs, uint32_t num, char **alleles)
{
  const char *first = seqs[0], *second = NULL;
  boolean more_than_two_alleles = false;
  size_t i;

  for(i = 1; i < num; i++) {
    if(strcmp(first, seqs[i]) != 0) {
      if(second == NULL) second = seqs[i];
      else if(strcmp(second, seqs[i]) != 0) {
        more_than_two_alleles = true;
        break;
      }
    }
  }

  // for(i = 0; i < num; i++) printf("%zu: %s\n", i, seqs[i]);

  assert(!more_than_two_alleles);

  // if(more_than_two_alleles)
  // {
  //   size_t msa_total_len = 0;
  //   multiple_seq_align(seqs, num, msa_wrking, &msa_total_len);

  //   size_t offset, msa_len = msa_total_len / num;

  //   alleles[0] = msa_wrking;
  //   for(i = 1, offset = msa_len; i < num; i++, offset += msa_len)
  //     alleles[i] = msa_wrking+offset;

  //   return msa_len;
  // }
  // else
  // {
    if(second == NULL) die("All alleles match");

    needleman_wunsch_align2(first, second, strlen(first), strlen(second),
                            &scoring, nw_aligner, alignment);

    for(i = 0; i < num; i++) {
      alleles[i] = strcmp(seqs[i], first) == 0 ? alignment->result_a
                                               : alignment->result_b;
    }

    return alignment->length;
  // }
}

static boolean is_allele_duplicate(vcf_entry_t *vcf, const StrBuf *allele)
{
  size_t i;
  if(strcmp(vcf->cols[VCFREF].buff, allele->buff) == 0) return true;
  for(i = 0; i < vcf->num_alts; i++)
    if(strcmp(vcf->alts[i].buff, allele->buff) == 0)
      return true;
  return false;
}

static void remove_duplicate_alleles(vcf_entry_t *vcf)
{
  size_t i = 0, num_alt_alleles = vcf->num_alts;
  StrBuf tmpbuf;
  vcf->num_alts = 0;

  while(i < num_alt_alleles)
  {
    if(is_allele_duplicate(vcf, &vcf->alts[i]))
    {
      SWAP(vcf->alts[i], vcf->alts[num_alt_alleles-1], tmpbuf);
      num_alt_alleles--;
    }
    else
    {
      i++;
    }
  }
  vcf->num_alts = i;
}

//
// Decompose into variants
//
static void parse_alignment(char **alleles, size_t num_alleles, size_t msa_len,
                            vcf_entry_t *outvcf, //size_t kmer_size,
                            const read_t *chr, size_t refpos)
{
  uint32_t allele_lens[MAX_ALLELES];

  // Keep track of which kmer we're on
  int start_idx[MAX_ALLELES];
  memset(start_idx, 0, num_alleles*sizeof(int));

  // Find where alleles differ
  int start, end = 0, j;
  size_t i;

  while((start = get_var_start(alleles, num_alleles, msa_len, end)) != -1)
  {
    // Update allele offsets
    for(i = 0; i < num_alleles; i++) {
      for(j = end; j < start; j++)
        if(alleles[i][j] != '-') start_idx[i]++;
    }

    end = get_var_end(alleles, num_alleles, msa_len, start);

    // REF position
    int pos = refpos + start_idx[0];

    strip_allele(&outvcf->cols[VCFREF], chr->seq.b + pos, end-start);
    boolean is_snp = (outvcf->cols[VCFREF].len == 1);
    allele_lens[0] = outvcf->cols[VCFREF].len;

    outvcf->num_alts = num_alleles - 1;

    for(i = 0; i < outvcf->num_alts; i++) {
      strip_allele(&outvcf->alts[i], alleles[i+1]+start, end-start);
      if(outvcf->alts[i].len != 1) is_snp = false;
      allele_lens[i+1] = outvcf->alts[i].len;
    }

    if(!is_snp)
    {
      char padding = pos > 0 ? chr->seq.b[pos-1] : 'N';
      strbuf_insert(&outvcf->cols[VCFREF], 0, &padding, 1);
      for(i = 0; i < outvcf->num_alts; i++) {
        strbuf_insert(&outvcf->alts[i], 0, &padding, 1);
      }
      pos--;
    }

    // Collapse down matching alleles
    remove_duplicate_alleles(outvcf);

    // printf("num_reduced_alleles: %zu; num_unused_alleles: %zu\n",
    //        num_reduced_alleles, num_unused_alleles);

    // Genotype
    /*
    uint32_t s, covgs[MAX_ALLELES];
    for(s = 0; s < num_samples; s++)
    {
      memset(covgs, 0, sizeof(uint32_t) * num_reduced_alleles);

      // Assign covg to allele
      int first_kmer, last_kmer;

      for(i = 0; i < outvcf->num_alts; i++)
      {
        first_kmer = start_idx[i+1] - kmer_size;
        last_kmer = start_idx[i+1] + allele_lens[i+1] - 1;

        first_kmer = MAX2(first_kmer, 0);
        first_kmer = MIN2(first_kmer, (signed)outvcf->covgs[s][i].len);
        last_kmer = MIN2(last_kmer, (signed)outvcf->covgs[s][i].len);

        covgs[allele_idx[i]]
          += supernode_read_starts(outvcf->covgs[s][i].arr + first_kmer,
                               last_kmer - first_kmer + 1);
      }

      StrBuf *genotype = &outvcf->cols[VCFSAMPLES+s];
      strbuf_reset(genotype);
      strbuf_sprintf(genotype, "%u", covgs[0]);
      for(i = 1; i < num_reduced_alleles; i++)
        strbuf_sprintf(genotype, ",%u", covgs[i]);
    }
    */

    // DEV: add pop filter here

    // Set ref position
    strbuf_reset(&outvcf->cols[VCFPOS]);
    strbuf_sprintf(&outvcf->cols[VCFPOS], "%i", pos+1);

    print_entry(outvcf, stdout);

    // Update offsets
    for(i = 0; i <= outvcf->num_alts; i++) {
      start_idx[i] += allele_lens[i];
    }
  }
}

static char isbam(const char *path, boolean bam)
{
  size_t len = strlen(path);
  return (len >= 4 && strcasecmp(path+len-4, bam ? ".bam" : ".sam") == 0);
}

static chrom_t* chrom_new()
{
  if(num_chroms == chrom_capacity) {
    chrom_capacity *= 2;
    chroms = realloc2(chroms, chrom_capacity * sizeof(chrom_t*));
  }
  seq_read_alloc(&chroms[num_chroms].r);
  return &chroms[num_chroms++];
}

static void chrom_free_last()
{
  num_chroms--;
  seq_read_dealloc(&chroms[num_chroms].r);
}

// Create chrom->read genome hash
static void load_ref_genome(char **paths, size_t num_files)
{
  seq_file_t *f;
  size_t i;
  khiter_t k;
  int hret;

  genome = kh_init(ghash);
  num_chroms = 0;
  chrom_capacity = 1024;
  chroms = malloc2(chrom_capacity * sizeof(chrom_t*));

  for(i = 0; i < num_files; i++)
  {
    if((f = seq_open(paths[i])) == NULL) die("Cannot open file %s\n", paths[i]);

    while(1)
    {
      chrom_t *chr = chrom_new();
      read_t *r = &(chr->r);

      if(seq_read(f,r) <= 0) { chrom_free_last(); break; }
      fprintf(stderr, "Chromosome [%s length:%zu]\n", r->name.b, r->seq.end);

      seq_read_to_uppercase(r);
      seq_read_truncate_name(r);

      k = kh_put(ghash, genome, r->name.b, &hret);
      if(hret == 0) {
        warn("duplicate chromosome (take first only): '%s'", r->name.b);
        chrom_free_last();
      }
      else {
        kh_value(genome, k) = chr;
      }
    }

    seq_close(f);
  }

  if(num_chroms == 0)
    die("No reference loaded");

  // fprintf(stderr, "Finished loading reference genome\n");
}

static void parse_header(gzFile gzvcf, StrBuf *line, CmdArgs *cmd,
                         char **refpaths, uint32_t num_refpaths)
{
  sample_indx = kh_init(samplehash);

  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));
  uint32_t i;
  boolean printed_info = false;

  char cwd[PATH_MAX + 1];

  samples_capacity = 16;
  sample_names = malloc2(samples_capacity * sizeof(char*));
  sample_total_seq = malloc2(samples_capacity * sizeof(uint64_t));

  while(strbuf_gzreadline_nonempty(line, gzvcf) > 0 && line->buff[0] == '#')
  {
    strbuf_chomp(line);
    char *str = line->buff;

    if(!strncasecmp(str, "##fileformat=", 13)) {
      printf("##fileformat=VCFv4.1\n");
    }
    else if(!strncasecmp(str, "##fileDate=", 11)) {
      printf("##fileDate=%s\n", datestr);
    }
    else if(!strncasecmp(str, "##reference=", 12)) {
      printf("##reference=file://%s", refpaths[0]);
      for(i = 1; i < num_refpaths; i++) printf(":%s", refpaths[i]);
      fputc('\n', stdout);
    }
    else if(!strncasecmp(str, "##phasing=", 10)) {
      printf("##phasing=partial\n");
      printf("##placeCmd=%s\n", cmd->cmdline);
      if(futil_get_current_dir(cwd) != NULL) printf("##placeCwd=%s\n", cwd);
      printf("##placeDate=%s\n", datestr);
    }
    else if(!strncasecmp(str, "##SAMPLE=", 9))
    {
      if(num_samples == samples_capacity)
      {
        samples_capacity *= 2;
        sample_names = realloc2(sample_names, samples_capacity * sizeof(char*));
        sample_total_seq = realloc2(sample_total_seq, samples_capacity * sizeof(uint64_t));
      }

      // strlen("##SAMPLE=<ID=") = 13
      char *sampleid = str+13, *sid_end, *totalseqstr, *seq_end;

      if(strncasecmp(str, "##SAMPLE=<ID=", 13) != 0)
        die("Unexpected VCF sample header line: %s", str);

      if((sid_end = strchr(sampleid, ',')) == NULL &&
         (sid_end = strchr(sampleid, '>')) == NULL)
      {
        die("Unexpected VCF sample header line: %s", str);
      }

      *sid_end = '\0';
      sample_names[num_samples] = strdup(sampleid);

      unsigned long totalseq = 0;

      if((totalseqstr = strstr(sid_end, "totalseqloaded")) != NULL)
      {
        totalseq = strtoul(totalseqstr+strlen("totalseqloaded"), &seq_end, 10);

        if(seq_end == totalseqstr+15)
          die("Unexpected VCF sample header line: %s", str);
      }

      sample_total_seq[num_samples] = totalseq;

      fprintf(stderr, "Sample [%s seq:%lu]\n", sample_names[num_samples], totalseq);

      khiter_t k;
      int hret;
      k = kh_put(samplehash, sample_indx, sample_names[num_samples], &hret);
      if(hret == 0) die("VCF multiple samples with the sample name");
      kh_value(sample_indx, k) = num_samples;

      num_samples++;
      *sid_end = ',';
      printf("%s\n", str);
    }
    else if(!strncasecmp(str, "##format=", 9) && !printed_info)
    {
      // print additional info tag headers
      // print additional filter headers
      printf(
"##INFO=<ID=BUBREV,Number=0,Type=Flag,Description=\"Bubble mapped to minus strand\">\n"
"##INFO=<ID=MQ,Number=0,Type=Flag,Description=\"MAPQ score of 5p flank mapping\">\n"
"##FILTER=<ID=MAPQ,Description=\"Mapped with MAPQ (<%u)\">\n"
"##FILTER=<ID=NOSAM,Description=\"No SAM entry\">\n"
"##FILTER=<ID=NOMAP,Description=\"Unmapped in SAM\">\n"
"##FILTER=<ID=FMAP,Description=\"3 prime flank did not align well\">\n"
"##FILTER=<ID=MAXBR,Description=\"Too many branches (>%u)\">\n"
"##FILTER=<ID=POPERR,Description=\"Pop filter classified as sequencing error\">\n"
"##FILTER=<ID=POPREP,Description=\"Pop filter classified as repeat\">\n"
"%s\n", min_mapq, max_branches, str);
      printed_info = true;
    }
    else if(!strncasecmp(str, "#CHROM", 6))
    {
      // printf("%s\n", str);
      size_t num_columns = count_char(str, '\t')+1;
      if(num_columns != VCFSAMPLES + num_samples)
        die("Incorrect number of columns in VCF file [num_samples: %u]", num_samples);
    }
    else printf("%s\n", str);

    strbuf_reset(line);
  }
}

static void parse_entry(vcf_entry_t *invcf, const bam1_t *bam,
                        vcf_entry_t *outvcf)
{
  // Copy details to output vcf
  vcf_entry_cpy(outvcf, invcf, num_samples);

  // Set chromosome field
  strbuf_set(&outvcf->cols[VCFCHROM], bam_header->target_name[bam->core.tid]);

  // look up chromosome
  khiter_t k = kh_get(ghash, genome, outvcf->cols[VCFCHROM].buff);
  if(k == kh_end(genome)) {
    print_entry(outvcf, stderr);
    die("Cannot find chrom [%s]", outvcf->cols[VCFCHROM].buff);
  }

  chrom_t *chrom = kh_value(genome, k);
  read_t *chr = &(chrom->r);

  // check orientation, offset
  // if rev-orient, revcmp alleles and flanks
  int cigar2rlen = bam_cigar2rlen(bam->core.n_cigar, bam_get_cigar(bam));

  // Get right flank and kmer_size
  size_t kmer_size = outvcf->lf->len-3;
  size_t endfl_missing = kmer_size - (outvcf->rf->len-3);

  if(bam_is_rev(bam)) vcf_entry_revcmp(outvcf);

  const StrBuf *lf = outvcf->lf, *rf = outvcf->rf;

  if(rf->len <= 3) rf = lf;

  strbuf_reset(endflank);
  if(bam_is_rev(bam)) {
    strbuf_copy(endflank, 0, lf->buff+3, lf->len-3);
    strbuf_append_strn(endflank, rf->buff+3, endfl_missing);
  } else {
    strbuf_copy(endflank, 0, rf->buff+3, rf->len-3);
    strbuf_append_strn(endflank, lf->buff+3, endfl_missing);
  }

  #ifdef DEBUG
    printf(" lf:%s rf:%s\n", lf->buff, rf->buff);
    printf(" endflank: %s\n", endflank->buff);
  #endif

  // end is index after last char
  int ref_start, ref_end;
  size_t longest_allele = vcf_entry_longest_allele(outvcf);

  if(bam_is_rev(bam)) {
    ref_start = bam->core.pos - longest_allele - kmer_size - 10;
    ref_end = bam->core.pos + kmer_size;
  } else {
    ref_start = bam->core.pos + cigar2rlen - kmer_size;
    ref_end = ref_start + kmer_size + longest_allele + kmer_size + 10;
  }

  if(ref_start < 0) ref_start = 0;
  if(ref_end > (signed)chr->seq.end) ref_end = chr->seq.end;

  char *ref_region = chr->seq.b+ref_start;
  int search_len = ref_end - ref_start;

  size_t i;
  size_t ref_trim_left = 0, ref_trim_right = 0;
  size_t bub_trim_left = 0, bub_trim_right = 0;
  char save_ref_base;

  // Attempt to find perfect match for kmer
  save_ref_base = ref_region[search_len];
  ref_region[search_len] = '\0';
  char *kmer_match = strstr(ref_region, endflank->buff);
  if(kmer_match == ref_region)
  {
    // If exact match at first kmer, look for a second match
    char *tmpsearch = strstr(ref_region+1, endflank->buff);
    if(tmpsearch != NULL) kmer_match = tmpsearch;
  }
  ref_region[search_len] = save_ref_base;

  #ifdef DEBUG
    printf(" ref_region: %.*s\n", search_len, ref_region);
  #endif

  if(kmer_match != NULL)
  {
    // Found exact match
    if(bam_is_rev(bam)) ref_trim_left = kmer_match - ref_region;
    else ref_trim_right = search_len - (kmer_match - ref_region + kmer_size);

    #ifdef DEBUG
      printf("sr:%.*s\nsr:", search_len, ref_region);
      uint32_t pos = kmer_match - ref_region;
      for(i = 0; i < pos; i++) fputc('-', stdout);
      fputs(endflank->buff, stdout);
      for(i = pos+kmer_size; i < (unsigned)search_len; i++) fputc('-', stdout);
      fputc('\n', stdout);
    #endif
  }
  else
  {
    // Look for approximate match
    needleman_wunsch_align2(ref_region, endflank->buff, search_len, kmer_size,
                            &scoring, nw_aligner, alignment);

    char *r1 = alignment->result_a, *r2 = alignment->result_b;

    #ifdef DEBUG
      printf("nw:%s\nnw:%s\n", r1, r2);
    #endif

    // --aa--cc-cge
    // aa--ccd-dcge

    // Find position of first match
    if(bam_is_rev(bam))
    {
      for(i = 0; i < alignment->length && r1[i] != r2[i]; i++) {
        ref_trim_left += (r1[i] != '-');
        bub_trim_left += (r2[i] != '-');
      }
    }
    else
    {
      for(i = alignment->length-1; r1[i] != r2[i]; i--) {
        ref_trim_right += (r1[i] != '-');
        bub_trim_right += (r2[i] != '-');
        if(i == 0) break;
      }
    }

    if(bub_trim_left + bub_trim_right > kmer_size / 2)
    {
      // flank doesn't map well
      vcf_entry_add_filter(outvcf, "FMAP");
      print_entry(outvcf, stdout);
      return;
    }
  }

  #ifdef DEBUG
    printf(" ref: %zu %zu; bub: %zu %zu\n",
           ref_trim_left, ref_trim_right,
           bub_trim_left, bub_trim_right);
  #endif

  // Check capacity of MSA array
  if(num_msa_alleles < outvcf->num_alts+1) {
    num_msa_alleles = ROUNDUP2POW(outvcf->num_alts+1);
    msa_alleles = realloc2(msa_alleles, num_msa_alleles * sizeof(char*));
  }

  char *ref_allele_str = chr->seq.b+ref_start+ref_trim_left;
  size_t ref_allele_len = search_len - ref_trim_left - ref_trim_right;
  save_ref_base = ref_allele_str[ref_allele_len];
  ref_allele_str[ref_allele_len] = '\0';

  size_t sum_len = 0;

  msa_alleles[0] = ref_allele_str;
  sum_len += ref_allele_len;

  // printf("lf: %s rf: %s [%zu:'%.*s']\n", lf->buff, rf->buff, rf->len-3-bub_trim_right,
  //        (int)(rf->len-3-bub_trim_right), rf->buff+3);

  for(i = 0; i < outvcf->num_alts; i++)
  {
    StrBuf *allele = &outvcf->alts[i];
    // printf("alleles: %s\n", allele->buff);
    strbuf_insert(allele, 0, lf->buff+3+bub_trim_left, lf->len-3-bub_trim_left);
    strbuf_append_strn(allele, rf->buff+3, rf->len-3-bub_trim_right);
    msa_alleles[i+1] = allele->buff;
    sum_len += allele->len;
  }

  if(sum_len+1 > msa_capacity) {
    msa_capacity = ROUNDUP2POW(sum_len+1);
    msa_wrking = realloc2(msa_wrking, msa_capacity * sizeof(char));
  }

  // Flanks not needed anymore
  strbuf_reset(outvcf->lf);
  strbuf_reset(outvcf->rf);

  size_t refpos = ref_allele_str - chr->seq.b;

  /*
  size_t msa_len = do_msa(msa_alleles, outvcf->num_alts+1, msa_alleles);

  #ifdef DEBUG
    for(i = 0; i < outvcf->num_alts+1; i++)
      printf(" msa: %.*s\n", (int)msa_len, msa_alleles[i]);
  #endif

  parse_alignment(msa_alleles, msa_len, outvcf, kmer_size, chr, refpos);
  */

  // for(i = 0; i <= outvcf->num_alts; i++)
  //   printf("%c) %s\n", (char)('a'+i), msa_alleles[i]);

  // print_entry(invcf, stdout);

  // Force pairwise alignment to ref
  // char *tmpalleles[2];
  // tmpalleles[0] = ref_allele_str;
  for(i = 1; i <= invcf->num_alts; i++) {
    if(strcmp(msa_alleles[i], ref_allele_str) != 0) {
      msa_alleles[0] = ref_allele_str;
      msa_alleles[1] = msa_alleles[i];
      size_t msa_len = do_msa(msa_alleles, 2, msa_alleles);
      // printf("X:%s\nY:%s\n", msa_alleles[0], msa_alleles[1]);
      parse_alignment(msa_alleles, 2, msa_len, outvcf, chr, refpos);
      // printf("//\n");
    }
  }

  // Re-instate removed ref base char
  ref_allele_str[ref_allele_len] = save_ref_base;
}

static uint32_t get_chrom_index(const char *chromname)
{
  khiter_t k = kh_get(ghash, genome, chromname);
  if(k == kh_end(genome))
    print_usage(usage, "No chrom with name [%s]", chromname);
  chrom_t *chrom = kh_value(genome, k);
  return chrom->index;
}

static uint32_t get_sample_index(const char *sample)
{
  khiter_t k = kh_get(samplehash, sample_indx, sample);
  if(k == kh_end(sample_indx))
    print_usage(usage, "No sample with name [%s]", sample);
  return kh_value(sample_indx, k);
}

static uint32_t check_ploidy(char *arg, char **sample, char **chrom)
{
  char *colon1, *colon2;
  uint32_t cpynum;

  if((colon1 = strchr(arg,':')) == NULL ||
     (colon2 = strchr(colon1+1,':')) == NULL ||
     colon1 == arg || colon2 == colon1+1 || colon2[1] == '\0' ||
     !parse_entire_uint(colon2+1, &cpynum))
  {
    print_usage(usage, "Invalid ploidy: %s", arg);
  }

  // Check sample name is valid
  *colon1 = '\0';
  get_sample_index(arg);
  *colon1 = ':';

  if(sample != NULL) {
    *colon1 = *colon2 = '\0';
    *sample = arg;
    *chrom = colon1+1;
  }

  return cpynum;
}

static void parse_ploidy(char *arg)
{
  char *samples, *chroms;
  uint32_t copy_num = check_ploidy(arg, &samples, &chroms);
  size_t i, end, chrindx, sindx, product = num_samples * num_chroms;

  if(strcmp(samples, "*") == 0) {
    if(strcmp(chroms, "*") == 0) {
      // All samples, all chromosomes
      for(i = 0; i < product; i++) ploidy[i] = copy_num;
    }
    else {
      char *chr = strtok(chroms, ",");
      do
      {
        // All samples, chromosome chrindx
        chrindx = get_chrom_index(chr);
        for(i = chrindx*num_samples, end = i+num_samples; i < end; i++)
          ploidy[i++] = copy_num;
      }
      while((chr = strtok(NULL, ",")) != NULL);
    }
  }
  else {
    char *sample = strtok(samples, ",");
    do
    {
      sindx = get_sample_index(sample);

      if(strcmp(chroms, "*") == 0)
      {
        // Sample sindx, all chromosomes
        for(i = sindx; i < product; i += num_samples) ploidy[i] = copy_num;
      }
      else {
        char *chr = strtok(chroms, ",");
        do
        {
          // Sample sindx, chromosome chrindx
          chrindx = get_chrom_index(chr);
          ploidy[(chrindx) * num_samples + (sindx)] = copy_num;
        }
        while((chr = strtok(NULL, ",")) != NULL);
      }
    } while((sample = strtok(NULL, ",")) != NULL);
  }
}

int ctx_place(CmdArgs *args)
{
  cmd_accept_options(args, "");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  // double x = gsl_sf_lnbeta(3,5);
  // printf("Result: %f\n", x);
  // exit(-1);

  // Use htslib
  // vcfFile *htsvcf = vcf_open(argv[1], "", 0);
  // bcf_hdr_t *header = vcf_hdr_read(htsvcf);
  // bcf1_t *entry = bcf_init1();
  // while(vcf_read1(htsvcf, header, entry) != -1) {
  //   printf("read pos:%i rlen:%i\n", (int)entry->pos, (int)entry->rlen);
  // }
  // bcf_destroy1(entry);
  // vcf_close(htsvcf);
  // exit(EXIT_FAILURE);

  int argi = 0;
  size_t i;

  // Read arguments
  while(argv[argi][0] == '-') {
    if(strcmp(argv[argi], "--minmapq") == 0) {
      if(argi + 1 == argc) print_usage(usage, NULL);
      if(!parse_entire_uint(argv[argi+1], &min_mapq))
        print_usage(usage, "Invalid --minmapq arg: %s", argv[argi+1]);
      argi++;
    } else if(strcmp(argv[argi], "--maxbr") == 0) {
      if(argi + 1 == argc) print_usage(usage, NULL);
      if(!parse_entire_uint(argv[argi+1], &max_branches))
        print_usage(usage, "Invalid --maxbr arg: %s", argv[argi+1]);
      if(max_branches > MAX_ALLELES)
        print_usage(usage, "--maxbr cannot be >%i", MAX_ALLELES);
      if(max_branches < 2)
        print_usage(usage, "--maxbr cannot be <2");
      argi++;
    } else if(strcmp(argv[argi], "--ploidy") == 0) {
      if(argi + 1 == argc) print_usage(usage, NULL);
      // check_ploidy(argv[argi+1], NULL, NULL);
      argi++;
    } else if(strcmp(argv[argi], "--genome") == 0) {
      if(argi + 1 == argc) print_usage(usage, NULL);
      if(!bases_to_integer(argv[argi+1], &genome_size)) {
        print_usage(usage, "Invalid <bases> arg (try 3.1G, 4.6Mbp etc): %s",
                    argv[argi+1]);
      }
      argi++;
    } else print_usage(usage, "Unknown argument: %s", argv[argi]);
    argi++;
  }

  if(argi+3 < argc) print_usage(usage, NULL);

  char *vcf_path = argv[argi++];
  char *sam_path = argv[argi++];

  char **ref_paths = argv + argi;
  size_t num_refpaths = argc - argi;

  // Check ploidy
  for(argi = 1; argv[argi][0] == '-'; argi++) {
    if(strcmp(argv[argi], "--minmapq") == 0) argi++;
    else if(strcmp(argv[argi], "--maxbr") == 0) argi++;
    else if(strcmp(argv[argi], "--ploidy") == 0) {
      check_ploidy(argv[argi+1], NULL, NULL);
      argi++;
    }
  }

  gzFile vcf = gzopen(vcf_path, "r");
  if(vcf == NULL) print_usage(usage, "Cannot open VCF %s", vcf_path);

  if(!isbam(sam_path, true) && !isbam(sam_path, false))
    print_usage(usage, "Mapped flanks is not .sam or .bam file: %s", sam_path);

  samFile *samfh = sam_open(sam_path, isbam(sam_path, true) ? "rb" : "rs", 0);
  if(samfh == NULL) die("Cannot open SAM/BAM %s", sam_path);

  // Load BAM header
  bam_header = sam_hdr_read(samfh);

  // Load VCF header
  StrBuf *line = strbuf_new();
  parse_header(vcf, line, args, ref_paths, num_refpaths);

  // Load reference genome
  load_ref_genome(ref_paths, num_refpaths);

  // Print remainder of VCF header
  for(i = 0; i < num_chroms; i++) {
    printf("##contig=<ID=%s,length=%zu>\n",
           chroms[i].r.name.b, chroms[i].r.seq.end);
  }
  printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  for(i = 0; i < num_samples; i++) {
    printf("\t%s", sample_names[i]);
  }
  printf("\n");

  // Now we have samples (from VCF header) and chroms (from reference genome)
  // Parse ploidy
  size_t product = num_samples * num_chroms;
  ploidy = malloc2(product * sizeof(uint32_t));

  // default to diploid
  for(i = 0; i < product; i++) ploidy[i] = 2;

  for(argi = 1; argv[argi][0] == '-'; argi++) {
    if(strcmp(argv[argi], "--minmapq") == 0) argi++;
    else if(strcmp(argv[argi], "--maxbr") == 0) argi++;
    else if(strcmp(argv[argi], "--ploidy") == 0) {
      parse_ploidy(argv[argi+1]);
      argi++;
    }
  }

  // print ploidy table to STDERR
  uint32_t indx = 0, j;
  fprintf(stderr, "ploidy:\n");
  for(j = 0; j < num_samples; j++) fprintf(stderr, "\t%s", sample_names[j]);
  fprintf(stderr, "\n");
  for(i = 0; i < num_chroms; i++) {
    fprintf(stderr, "[%s]", chroms[i].r.name.b);
    for(j = 0; j < num_samples; j++, indx++) {
      fprintf(stderr, "\t%u", ploidy[indx]);
    }
    fprintf(stderr, "\n");
  }

  // Setup for loading VCF lines
  vcf_entry_t invcf, outvcf;
  vcf_entry_alloc(&invcf, num_samples);
  vcf_entry_alloc(&outvcf, num_samples);

  // Setup pairwise aligner
  nw_aligner = needleman_wunsch_new();
  alignment = alignment_create(1024);
  scoring_init(&scoring, 1, -2, -4, -1, true, true, 0, 0, 0, 0);

  // Setup multiple sequence alignment
  msa_wrking = malloc2(msa_capacity * sizeof(char));
  msa_alleles = malloc2(num_msa_alleles * sizeof(char*));
  endflank = strbuf_new();

  bam1_t *bam = bam_init1();
  boolean read_sam = (sam_read1(samfh, bam_header, bam) >= 0);
  boolean read_vcf = (line->len > 0);

  while(read_sam && read_vcf)
  {
    vcf_entry_parse(line, &invcf, num_samples);

    if(strcmp(invcf.cols[VCFCHROM].buff, "un") != 0 ||
       strcmp(invcf.cols[VCFPOS].buff, "1") != 0 ||
       strcmp(invcf.cols[VCFREF].buff, "N") != 0)
    {
      die("Unexpected VCF line: %s", line->buff);
    }

    // Get bam query name
    char *bname = bam_get_qname(bam);

    if(strcmp(invcf.cols[VCFID].buff, bname) != 0)
    {
      // BAM entry and VCF name do not match - skip
      vcf_entry_add_filter(&invcf, "NOSAM");
      print_entry(&invcf, stdout);

      // Read next vcf entry only (not sam)
      strbuf_reset(line);
      read_vcf = strbuf_gzreadline_nonempty(line, vcf) > 0;
    }
    else
    {
      int unmapped = bam->core.flag & BAM_FUNMAP;
      int low_mapq = bam->core.qual < min_mapq;
      int toomanybr = invcf.num_alts > max_branches;

      if(!unmapped) {
        info_tag_add(&invcf, "MQ=%i", (int)bam->core.qual);
      }

      if(unmapped || low_mapq || toomanybr)
      {
        // Entry failed a filter
        if(unmapped) vcf_entry_add_filter(&invcf, "NOMAP");
        if(low_mapq) vcf_entry_add_filter(&invcf, "MAPQ");
        if(toomanybr) vcf_entry_add_filter(&invcf, "MAXBR");

        print_entry(&invcf, stdout);
      }
      else if(invcf.lf == NULL || invcf.rf == NULL) {
        print_entry(&invcf, stderr);
        die("Missing LF/RF tags");
      }
      else {
        parse_entry(&invcf, bam, &outvcf);
      }

      // Read next sam and vcf
      read_sam = (sam_read1(samfh, bam_header, bam) >= 0);
      strbuf_reset(line);
      read_vcf = strbuf_gzreadline_nonempty(line, vcf) > 0;
    }
  }

  if(read_sam) warn("Excess sam entries: %s", bam_get_qname(bam));
  if(read_vcf) warn("Excess vcf entries: %s", line->buff);

  free(msa_wrking);
  free(msa_alleles);

  strbuf_free(endflank);

  kh_destroy(ghash, genome);
  kh_destroy(samplehash, sample_indx);

  for(i = 0; i < num_samples; i++) free(sample_names[i]);
  free(sample_names);
  free(sample_total_seq);

  while(num_chroms > 0) chrom_free_last();
  free(chroms);

  alignment_free(alignment);
  needleman_wunsch_free(nw_aligner);

  sam_close(samfh);
  free(bam_header);
  free(bam);

  strbuf_free(line);
  vcf_entry_dealloc(&invcf, num_samples);
  vcf_entry_dealloc(&outvcf, num_samples);

  return EXIT_SUCCESS;
}
