#include "global.h"
#include <time.h>

// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_cdf.h> // beta distribution
// #include <gsl/gsl_sf_gamma.h> // gamma and beta functions

// #include "vcf.h" // not using htslib
#include "sam.h"
#include "khash.h"
#include "seq_file.h"
#include "string_buffer.h"
#include "needleman_wunsch.h"

#include "tools.h"
#include "util.h"
#include "file_util.h"
#include "vcf_parsing.h"
#include "binary_kmer.h"
#include "supernode.h"

// #define CTXVERBOSE 1

int nwmatch = 1, nwmismatch = -2, nwgapopen = -4, nwgapextend = -1;

const char place_usage[] =
"usage: "CMD" place [options] <calls.vcf> <calls.sam> <ref1.fa ...>\n"
"  Align calls to a reference genome.\n"
"  Options:\n"
"    --minmapq <mapq>   Flank must map with MAPQ >= <mapq> [default: 30]\n"
"    --out <out.vcf>    Output file [default: STDOUT]\n"
"\n"
"  Alignment:\n"
"    --match <m>      [default:  1]\n"
"    --mismatch <m>   [default: -2]\n"
"    --gapopen <m>    [default: -4]\n"
"    --gapextend <m>  [default: -1]\n";

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
size_t num_samples = 0, samples_capacity;
char **sample_names;
uint64_t *sample_total_seq;

// nw alignment
nw_aligner_t *nw_aligner;
alignment_t *alignment;
scoring_t *nw_scoring_flank, *nw_scoring_allele;
size_t num_nw_flank = 0, num_nw_allele = 0;

// Temporary memory
StrBuf endflank;

// Filtering parameters
size_t min_mapq = 30;

// VCF printing
size_t num_variants_printed = 0;

static void print_entry(vcf_entry_t *vcfentry, FILE *out)
{
  // Set var num
  strbuf_reset(&vcfentry->cols[VCFID]);
  strbuf_sprintf(&vcfentry->cols[VCFID], "var%zu", num_variants_printed++);
  vcf_entry_print(vcfentry, out, num_samples);
}

// Returns msa_len if no more variants
static size_t get_var_start(char **alleles, size_t num_alleles,
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

static size_t get_var_end(char **alleles, size_t num_alleles,
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

static void strip_allele(StrBuf *allele, const char* input, size_t len)
{
  strbuf_reset(allele);
  size_t i;
  for(i = 0; i < len; i++) {
    if(input[i] != '-') strbuf_append_char(allele, input[i]);
  }
}

static bool is_allele_duplicate(const vcf_entry_t *vcf, const StrBuf *allele)
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

  while(i < num_alt_alleles) {
    if(is_allele_duplicate(vcf, &vcf->alts[i])) {
      SWAP(vcf->alts[i], vcf->alts[num_alt_alleles-1], tmpbuf);
      num_alt_alleles--;
    }
    else i++;
  }
  vcf->num_alts = i;
}

//
// Decompose into variants
//
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

// Returns 0 if failed to find 3p flank, 1 otherwise
static int parse_entry(const vcf_entry_t *invcf, const bam1_t *bam,
                       vcf_entry_t *outvcf, FILE *fout)
{
  if(bam_is_rev(bam)) ctx_assert(invcf->lf->len <= invcf->rf->len);
  else ctx_assert(invcf->lf->len >= invcf->rf->len);

  // Copy details to output vcf and set chromosome field
  const char *chrname = bam_header->target_name[bam->core.tid];
  vcf_entry_cpy(outvcf, invcf, num_samples);
  strbuf_set(&outvcf->cols[VCFCHROM], chrname);

  // look up chromosome
  khiter_t k = kh_get(ghash, genome, chrname);
  if(k == kh_end(genome)) {
    print_entry(outvcf, stderr);
    die("Cannot find chrom [%s]", chrname);
  }

  const chrom_t *chrom = kh_value(genome, k);
  const read_t *chr = &(chrom->r);

  // check orientation, offset
  // if rev-orient, revcmp alleles and flanks
  int cigar2rlen = bam_cigar2rlen(bam->core.n_cigar, bam_get_cigar(bam));

  // Get right flank and kmer_size
  // Both left and right flanks are prefixed with LF= and RF= respectively
  // hence lots of +3 etc
  const size_t FPREFIX = 3; // strlen("LF=")
  const size_t kmer_size = MAX2(invcf->lf->len,invcf->rf->len)-FPREFIX;
  ctx_assert(kmer_size >= 3);

  size_t endfl_missing;

  const StrBuf *lf, *rf;
  lf = (invcf->lf->len <= FPREFIX ? invcf->rf : invcf->lf);
  rf = (invcf->rf->len <= FPREFIX ? invcf->lf : invcf->rf);

  strbuf_reset(&endflank);

  if(bam_is_rev(bam)) {
    endfl_missing = kmer_size - (lf->len-FPREFIX);
    strbuf_append_strn(&endflank, lf->buff+FPREFIX, lf->len-FPREFIX);
    strbuf_append_strn(&endflank, rf->buff+FPREFIX, endfl_missing);
  } else {
    endfl_missing = kmer_size - (rf->len-FPREFIX);
    strbuf_append_strn(&endflank, lf->buff+lf->len-endfl_missing, endfl_missing);
    strbuf_append_strn(&endflank, rf->buff+FPREFIX, rf->len-FPREFIX);
  }

  #ifdef CTXVERBOSE
    printf(" %s\n", bam_is_rev(bam) ? "reverse" : "forward");
    printf(" lf:%s rf:%s\n", lf->buff, rf->buff);
    printf(" endflank: %s\n", endflank.buff);
  #endif

  // Choose a region of the ref to search for the end flank
  // end is index after last char
  long search_start, search_end;
  size_t longest_allele = vcf_entry_longest_allele(invcf);

  if(bam_is_rev(bam)) {
    search_start = (long)bam->core.pos - (long)(longest_allele + kmer_size + 10);
    search_end = (long)bam->core.pos + (long)kmer_size;
  } else {
    search_start = (long)(bam->core.pos + cigar2rlen) - (long)kmer_size;
    search_end = search_start + (long)(kmer_size + longest_allele + kmer_size + 10);
  }

  #ifdef CTXVERBOSE
    printf(" search start:%li end:%li\n", search_start, search_end);
  #endif

  if(search_start < 0) search_start = 0;
  if(search_end > (signed)chr->seq.end) search_end = (long)chr->seq.end;

  char *search_region = chr->seq.b+search_start;
  size_t search_len = (size_t)(search_end - search_start);

  size_t i;
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

  #ifdef CTXVERBOSE
    printf(" search_region: %.*s\n", (int)search_len, search_region);
  #endif

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
    num_nw_flank++;

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

    if(matches < kmer_size / 2)
    {
      // flank doesn't map well
      vcf_entry_add_filter(outvcf, "FMAP");
      print_entry(outvcf, fout);
      return 0;
    }
  }

  #ifdef CTXVERBOSE
    printf(" ref: %zu %zu; bub: %zu %zu\n",
           search_trim_left, search_trim_right,
           flank_trim_left, flank_trim_right);
  #endif

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

  // Force pairwise alignment to ref
  // status("reflen: %zu\n", reflen);

  for(i = 0; i < invcf->num_alts; i++)
  {
    if(reflen != invcf->alts[i].len ||
       strcmp(ref_allele_str, invcf->alts[i].buff) != 0)
    {
      needleman_wunsch_align2(ref_allele_str, invcf->alts[i].buff,
                              reflen, invcf->alts[i].len,
                              nw_scoring_allele, nw_aligner, alignment);
      num_nw_allele++;

      #ifdef CTXVERBOSE
        printf("X:%s\nY:%s\n", alignment->result_a, alignment->result_b);
      #endif

      char *alignments[2] = {alignment->result_a, alignment->result_b};
      parse_alignment(alignments, 2, alignment->length, chr, refpos, outvcf, fout);
    }
  }

  return 1;
}

static bool isbam(const char *path, bool bam)
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
  if(seq_read_alloc(&chroms[num_chroms].r) == NULL) die("Out of memory");
  return &chroms[num_chroms++];
}

static void chrom_free_last()
{
  seq_read_dealloc(&chroms[--num_chroms].r);
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
      status("Chromosome [%s length:%zu]", r->name.b, r->seq.end);

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

  status("Finished loading reference genome");
}

static void parse_header(gzFile gzvcf, StrBuf *line, CmdArgs *cmd,
                         char *const* refpaths, size_t num_refpaths, FILE *fout)
{
  sample_indx = kh_init(samplehash);

  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));
  size_t i;
  bool printed_info = false;

  char cwd[PATH_MAX + 1];

  samples_capacity = 16;
  sample_names = malloc2(samples_capacity * sizeof(char*));
  sample_total_seq = malloc2(samples_capacity * sizeof(uint64_t));

  while(strbuf_gzreadline_nonempty(line, gzvcf) > 0 && line->buff[0] == '#')
  {
    strbuf_chomp(line);
    char *str = line->buff;

    if(!strncasecmp(str, "##fileformat=", 13)) {
      fprintf(fout, "##fileformat=VCFv4.1\n");
    }
    else if(!strncasecmp(str, "##fileDate=", 11)) {
      fprintf(fout, "##fileDate=%s\n", datestr);
    }
    else if(!strncasecmp(str, "##reference=", 12)) {
      fprintf(fout, "##reference=file://%s", refpaths[0]);
      for(i = 1; i < num_refpaths; i++) {
        fputc(':', fout); fputs(refpaths[i], fout);
      }
      fputc('\n', fout);
    }
    else if(!strncasecmp(str, "##phasing=", 10)) {
      fprintf(fout, "##phasing=partial\n");
      fprintf(fout, "##placeCmd=%s\n", cmd->cmdline);
      if(futil_get_current_dir(cwd) != NULL) fprintf(fout, "##placeCwd=%s\n", cwd);
      fprintf(fout, "##placeDate=%s\n", datestr);
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

      status("Sample [%s seq:%lu]\n", sample_names[num_samples], totalseq);

      khiter_t k;
      int hret;
      k = kh_put(samplehash, sample_indx, sample_names[num_samples], &hret);
      if(hret == 0) die("VCF multiple samples with the sample name");
      kh_value(sample_indx, k) = (uint32_t)num_samples;

      num_samples++;
      *sid_end = ',';
      fputs(str, fout);
      fputc('\n', fout);
    }
    else if(!strncasecmp(str, "##format=", 9) && !printed_info)
    {
      // print additional info tag headers
      // print additional filter headers
      fprintf(fout,
"##INFO=<ID=BUBREV,Number=0,Type=Flag,Description=\"Bubble mapped to minus strand\">\n"
"##INFO=<ID=MQ,Number=0,Type=Flag,Description=\"MAPQ score of 5p flank mapping\">\n"
"##FILTER=<ID=MAPQ,Description=\"Mapped with MAPQ (<%zu)\">\n"
"##FILTER=<ID=NOSAM,Description=\"No SAM entry\">\n"
"##FILTER=<ID=NOMAP,Description=\"Unmapped in SAM\">\n"
"##FILTER=<ID=FMAP,Description=\"3 prime flank did not align with >kmer/2 base matches\">\n"
"##FILTER=<ID=POPERR,Description=\"Pop filter classified as sequencing error\">\n"
"##FILTER=<ID=POPREP,Description=\"Pop filter classified as repeat\">\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
// "##FORMAT=<ID=DP,Number=.,Type=Integer,Description=\"Read Depth\">\n" // also Number=R
"%s\n", min_mapq, str);
      printed_info = true;
    }
    else if(!strncasecmp(str, "#CHROM", 6))
    {
      // fprintf(fout, "%s\n", str);
      size_t num_columns = count_char(str, '\t')+1;
      if(num_columns != VCFSAMPLES + num_samples)
        die("Incorrect number of columns in VCF file [num_samples: %zu]", num_samples);
    }
    else { fputs(str, fout); fputc('\n', fout); }

    strbuf_reset(line);
  }
}

int ctx_place(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Have already checked we have at least 3 arguments

  // double x = gsl_sf_lnbeta(3,5);
  // printf("Result: %f\n", x);
  // exit(-1);

  // Use htslib
  // vcfFile *htsvcf = vcf_open(argv[1], "r");
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

  const char *nwargs[4] = {"--match", "--mismatch", "--gapopen", "--gapextend"};
  int *nwargptrs[4] = {&nwmatch, &nwmismatch, &nwgapopen, &nwgapextend};

  // Read arguments
  while(argv[argi][0] == '-') {
    if(strcmp(argv[argi], "--minmapq") == 0) {
      if(argi + 1 == argc) cmd_print_usage(NULL);
      if(!parse_entire_size(argv[argi+1], &min_mapq))
        cmd_print_usage("Invalid --minmapq arg: %s", argv[argi+1]);
      argi++;
    }
    else
    {
      // Parse an alignment parameter
      for(i = 0; i < 4 && strcmp(argv[argi], nwargs[i]) != 0; i++);
      if(i == 4) cmd_print_usage("Unknown argument: %s", argv[argi]);
      if(argi + 1 == argc) cmd_print_usage(NULL);
      if(!parse_entire_int(argv[argi+1], nwargptrs[i]))
        cmd_print_usage("Invalid %s arg: %s", argv[argi], argv[argi+1]);
      argi++;
    }
    argi++;
  }

  if(argi+3 < argc) cmd_print_usage(NULL);

  // Check alignment args
  if(nwmatch < MAX3(nwmismatch, nwgapopen, nwgapextend)) {
    cmd_print_usage("Alignment match should be greater than "
                       "mismatch, gap open and extend");
  }

  char *vcf_path = argv[argi++];
  char *sam_path = argv[argi++];

  char **ref_paths = argv + argi;
  size_t num_refpaths = (size_t)(argc - argi);

  gzFile vcf = gzopen(vcf_path, "r");
  if(vcf == NULL) cmd_print_usage("Cannot open VCF %s", vcf_path);

  if(!isbam(sam_path, true) && !isbam(sam_path, false))
    cmd_print_usage("Mapped flanks is not .sam or .bam file: %s", sam_path);

  samFile *samfh = sam_open(sam_path, isbam(sam_path, true) ? "rb" : "rs");
  if(samfh == NULL) die("Cannot open SAM/BAM %s", sam_path);

  FILE *fout = stdout;
  if(args->output_file_set) {
    fout = fopen(args->output_file, "w");
    if(fout == NULL) cmd_print_usage("Cannot open file: %s", args->output_file);
  }

  // Load BAM header
  bam_header = sam_hdr_read(samfh);

  // Load VCF header
  StrBuf *line = strbuf_new();
  parse_header(vcf, line, args, ref_paths, num_refpaths, fout);

  // Load reference genome
  load_ref_genome(ref_paths, num_refpaths);

  // Print remainder of VCF header
  for(i = 0; i < num_chroms; i++) {
    fprintf(fout, "##contig=<ID=%s,length=%zu>\n",
            chroms[i].r.name.b, chroms[i].r.seq.end);
  }
  fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", fout);
  for(i = 0; i < num_samples; i++) {
    fputc('\t', fout); fputs(sample_names[i], fout);
  }
  fputc('\n', fout);

  // Setup for loading VCF lines
  vcf_entry_t invcf, outvcf;
  vcf_entry_alloc(&invcf, num_samples);
  vcf_entry_alloc(&outvcf, num_samples);

  // Setup pairwise aligner
  nw_aligner = needleman_wunsch_new();
  alignment = alignment_create(1024);
  nw_scoring_flank = malloc2(sizeof(scoring_t));
  nw_scoring_allele = malloc2(sizeof(scoring_t));
  scoring_init(nw_scoring_flank, nwmatch, nwmismatch, nwgapopen, nwgapextend,
               true, true, 0, 0, 0, 0);
  scoring_init(nw_scoring_allele, nwmatch, nwmismatch, nwgapopen, nwgapextend,
               false, false, 0, 0, 0, 0);

  status("Alignment match:%i mismatch:%i gapopen:%i gapextend:%i",
         nwmatch, nwmismatch, nwgapopen, nwgapextend);

  strbuf_alloc(&endflank, 1024);

  bam1_t *bam = bam_init1();
  bool read_sam = (sam_read1(samfh, bam_header, bam) >= 0);
  bool read_vcf = (line->len > 0);

  // Filter statistics
  size_t num_missing_sam = 0, num_unmapped = 0, num_low_mapq = 0;
  size_t num_3p_not_found = 0, num_passed = 0, num_bubbles = 0;

  while(read_sam && read_vcf)
  {
    vcf_entry_parse(line, &invcf, num_samples);
    num_bubbles++;

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
      print_entry(&invcf, fout);
      num_missing_sam++;

      // Read next vcf entry only (not sam)
      strbuf_reset(line);
      read_vcf = strbuf_gzreadline_nonempty(line, vcf) > 0;
    }
    else
    {
      int unmapped = bam->core.flag & BAM_FUNMAP;
      int low_mapq = bam->core.qual < min_mapq;

      if(!unmapped) vcf_info_tag_add(&invcf, "MQ=%i", (int)bam->core.qual);

      if(unmapped || low_mapq)
      {
        // Entry failed a filter
        if(unmapped) { vcf_entry_add_filter(&invcf, "NOMAP"); num_unmapped++; }
        if(low_mapq) { vcf_entry_add_filter(&invcf, "MAPQ"); num_low_mapq++; }

        print_entry(&invcf, fout);
      }
      else if(invcf.lf == NULL || invcf.rf == NULL) {
        print_entry(&invcf, stderr);
        die("Missing LF/RF tags");
      }
      else {
        if(bam_is_rev(bam)) vcf_entry_revcmp(&invcf);
        if(parse_entry(&invcf, bam, &outvcf, fout)) num_passed++;
        else num_3p_not_found++;
      }

      // Read next sam and vcf
      read_sam = (sam_read1(samfh, bam_header, bam) >= 0);
      strbuf_reset(line);
      read_vcf = strbuf_gzreadline_nonempty(line, vcf) > 0;
    }
  }

  if(read_sam) warn("Excess sam entries: %s", bam_get_qname(bam));
  if(read_vcf) warn("Excess vcf entries: %s", line->buff);

  char missing_sam_str[100], unmapped_str[100], minmapq_str[100];
  char flank3p_not_found_str[100];
  char passed_str[100], total_bubbles_str[100], nvariants_str[100];

  ulong_to_str(num_missing_sam, missing_sam_str);
  ulong_to_str(num_unmapped, unmapped_str);
  ulong_to_str(num_low_mapq, minmapq_str);
  ulong_to_str(num_3p_not_found, flank3p_not_found_str);
  ulong_to_str(num_passed, passed_str);
  ulong_to_str(num_bubbles, total_bubbles_str);
  ulong_to_str(num_variants_printed, nvariants_str);

  double pass_rate = 100.0 * (double)num_passed / num_bubbles;

  char nw_flank_str[100], nw_allele_str[100];
  ulong_to_str(num_nw_flank, nw_flank_str);
  ulong_to_str(num_nw_allele, nw_allele_str);

  status("Used %s NW for flanks, %s NW for alleles", nw_flank_str, nw_allele_str);
  status("Bubbles missing SAM entry: %s", missing_sam_str);
  status("Bubbles unmapped: %s", unmapped_str);
  status("Bubbles MAPQ < %zu: %s", min_mapq, minmapq_str);
  status("Bubbles 3p flank not found: %s", flank3p_not_found_str);
  status("Passed bubbles: %s / %s [%.2f%%]", passed_str, total_bubbles_str, pass_rate);
  status("Variants printed: %s", nvariants_str);
  status("Variants written to: %s",
         args->output_file_set ? args->output_file : "STDOUT");

  if(args->output_file_set) fclose(fout);

  strbuf_dealloc(&endflank);

  kh_destroy(ghash, genome);
  kh_destroy(samplehash, sample_indx);

  for(i = 0; i < num_samples; i++) free(sample_names[i]);
  free(sample_names);
  free(sample_total_seq);

  while(num_chroms > 0) chrom_free_last();
  free(chroms);

  free(nw_scoring_flank);
  free(nw_scoring_allele);
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
