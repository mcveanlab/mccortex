#include "global.h"
#include "seq_file.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "vcf_parsing.h"
#include "binary_kmer.h"
#include "supernode.h"
#include "seq_reader.h"

// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_cdf.h> // beta distribution
// #include <gsl/gsl_sf_gamma.h> // gamma and beta functions

// #include "vcf.h" // not using htslib
#include "sam.h"
#include "seq-align/src/needleman_wunsch.h"

#include <time.h>

// #define CTXVERBOSE 1

const char place_usage[] =
"usage: "CMD" place [options] <calls.vcf> <calls.sam> <ref1.fa ...>\n"
"\n"
"  Align calls to a reference genome.\n"
"\n"
"  -h, --help             This help message\n"
"  -o, --out <out.vcf>    Output file [default: STDOUT]\n"
"  -Q, --minmapq <mapq>   Flank must map with MAPQ >= <mapq> [default: 30]\n"
"\n"
"  Alignment scoring:\n"
"  -m, --match <m>      [default:  1]\n"
"  -M, --mismatch <m>   [default: -2]\n"
"  -g, --gapopen <m>    [default: -4]\n"
"  -G, --gapextend <m>  [default: -1]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
// command specific
  {"minmapq",      required_argument, NULL, 'Q'},
// alignment
  {"match",        required_argument, NULL, 'm'},
  {"mismatch",     required_argument, NULL, 'M'},
  {"gap-open",     required_argument, NULL, 'g'},
  {"gap-extend",   required_argument, NULL, 'G'},
  {NULL, 0, NULL, 0}
};

const char *out_path = NULL;

// Reference genome
static khash_t(ChromHash) *genome;
static ReadBuffer chroms;

// Flank mapping
static bam_hdr_t *bam_header;

// VCF header info
KHASH_MAP_INIT_STR(samplehash, uint32_t)
static khash_t(samplehash) *sample_indx;
static size_t num_samples = 0, samples_capacity;
static char **sample_names;
static uint64_t *sample_total_seq;

// nw alignment
static int nwmatch = 1, nwmismatch = -2, nwgapopen = -4, nwgapextend = -1;
static nw_aligner_t *nw_aligner;
static alignment_t *alignment;
static scoring_t *nw_scoring_flank, *nw_scoring_allele;
static size_t num_nw_flank = 0, num_nw_allele = 0;

// Temporary memory
static StrBuf endflank;

// Filtering parameters
#define DEFAULT_MIN_MAPQ 30
static size_t min_mapq = SIZE_MAX;

// VCF printing
static size_t num_variants_printed = 0;

static void print_entry(vcf_entry_t *vcfentry, FILE *out)
{
  // Set var num
  strbuf_reset(&vcfentry->cols[VCFID]);
  strbuf_sprintf(&vcfentry->cols[VCFID], "var%zu", num_variants_printed++);
  vcf_entry_print(vcfentry, out, 0);
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
  vcf->num_alts = 0;

  while(i < num_alt_alleles) {
    if(is_allele_duplicate(vcf, &vcf->alts[i])) {
      SWAP(vcf->alts[i], vcf->alts[num_alt_alleles-1]);
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
  khiter_t k = kh_get(ChromHash, genome, chrname);
  if(k == kh_end(genome)) {
    print_entry(outvcf, stderr);
    die("Cannot find chrom [%s]", chrname);
  }

  const read_t *chr = kh_value(genome, k);

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

static void parse_header(gzFile gzvcf, StrBuf *line,
                         char **refpaths, size_t num_ref_paths,
                         FILE *fout)
{
  sample_indx = kh_init(samplehash);

  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));
  size_t i;

  samples_capacity = 16;
  sample_names = ctx_malloc(samples_capacity * sizeof(char*));
  sample_total_seq = ctx_malloc(samples_capacity * sizeof(uint64_t));

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
      for(i = 1; i < num_ref_paths; i++) {
        fputc(':', fout); fputs(refpaths[i], fout);
      }
      fputc('\n', fout);
    }
    else if(!strncasecmp(str, "##phasing=", 10)) {
      fprintf(fout, "##phasing=partial\n");
      fprintf(fout, "##placeCmd=%s\n", cmd_get_cmdline());
      fprintf(fout, "##placeCwd=%s\n", cmd_get_cwd());
      fprintf(fout, "##placeDate=%s\n", datestr);
    }
    else if(!strncasecmp(str, "##colour=", 9))
    {
      if(num_samples == samples_capacity)
      {
        samples_capacity *= 2;
        sample_names = ctx_realloc(sample_names, samples_capacity * sizeof(char*));
        sample_total_seq = ctx_realloc(sample_total_seq, samples_capacity * sizeof(uint64_t));
      }

      // strlen("##colour=<ID=") = 13
      char *sampleid = str+13, *sid_end, *totalseqstr, *seq_end;

      if(strncasecmp(str, "##colour=<ID=", 13) != 0)
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
    else if(!strncasecmp(str, "#CHROM", 6))
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
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", min_mapq);

      size_t num_columns = string_count_char(str, '\t')+1;
      if(num_columns != VCFSAMPLES)
        die("Incorrect number of columns in VCF file [expect: %i]", VCFSAMPLES);
    }
    else { fputs(str, fout); fputc('\n', fout); }

    strbuf_reset(line);
  }
}

int ctx_place(int argc, char **argv)
{
  // hide unused function warnings
  (void)kh_clear_ChromHash;
  (void)kh_del_ChromHash;
  (void)kh_clear_samplehash;
  (void)kh_get_samplehash;
  (void)kh_del_samplehash;

  size_t i;

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

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': cmd_check(out_path != NULL, cmd); out_path = optarg; break;
      case 'Q': cmd_check(min_mapq != SIZE_MAX,cmd); min_mapq = cmd_uint32(cmd, optarg); break;
      case 'g': nwmatch = cmd_int32(cmd, optarg); break;
      case 'G': nwmismatch = cmd_int32(cmd, optarg); break;
      case 'm': nwgapopen = cmd_int32(cmd, optarg); break;
      case 'M': nwgapextend = cmd_int32(cmd, optarg); break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" place -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(min_mapq == SIZE_MAX) min_mapq = DEFAULT_MIN_MAPQ;

  if(optind+3 < argc) cmd_print_usage(NULL);

  // Check alignment args
  if(nwmatch < MAX3(nwmismatch, nwgapopen, nwgapextend)) {
    cmd_print_usage("Alignment match should be greater than "
                       "mismatch, gap open and extend");
  }

  char *vcf_path = argv[optind++];
  char *sam_path = argv[optind++];

  char **ref_paths = argv + optind;
  size_t num_ref_paths = (size_t)(argc - optind);

  gzFile vcf = gzopen(vcf_path, "r");
  if(vcf == NULL) cmd_print_usage("Cannot open VCF %s", vcf_path);

  if(!futil_path_has_extension(sam_path, ".bam") &&
     !futil_path_has_extension(sam_path, ".sam"))
  {
    cmd_print_usage("Mapped flanks is not .sam or .bam file: %s", sam_path);
  }

  bool isbam = futil_path_has_extension(sam_path, ".bam");

  samFile *samfh = sam_open(sam_path, isbam ? "rb" : "rs");
  if(samfh == NULL) die("Cannot open SAM/BAM %s", sam_path);

  FILE *fout = futil_open_output(out_path);

  // Load BAM header
  bam_header = sam_hdr_read(samfh);

  // Load VCF header
  StrBuf *line = strbuf_new();
  parse_header(vcf, line, ref_paths, num_ref_paths, fout);

  // Load reference genome
  readbuf_alloc(&chroms, 1024);
  genome = kh_init(ChromHash);
  seq_reader_load_ref_genome(ref_paths, num_ref_paths, &chroms, genome);

  // Print remainder of VCF header
  for(i = 0; i < chroms.len; i++) {
    fprintf(fout, "##contig=<ID=%s,length=%zu>\n",
            chroms.data[i].name.b, chroms.data[i].seq.end);
  }
  fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n", fout);

  // Setup for loading VCF lines
  vcf_entry_t invcf, outvcf;
  vcf_entry_alloc(&invcf, num_samples);
  vcf_entry_alloc(&outvcf, num_samples);

  // Setup pairwise aligner
  nw_aligner = needleman_wunsch_new();
  alignment = alignment_create(1024);
  nw_scoring_flank = ctx_malloc(sizeof(scoring_t));
  nw_scoring_allele = ctx_malloc(sizeof(scoring_t));
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
    vcf_entry_parse(line, &invcf, 0);
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

  double pass_rate = num_bubbles ? 100.0 * (double)num_passed / num_bubbles : 0;

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
  status("Variants written to: %s", futil_outpath_str(out_path));

  fclose(fout);

  strbuf_dealloc(&endflank);

  kh_destroy(ChromHash, genome);
  kh_destroy(samplehash, sample_indx);

  for(i = 0; i < num_samples; i++) ctx_free(sample_names[i]);
  ctx_free(sample_names);
  ctx_free(sample_total_seq);

  for(i = 0; i < chroms.len; i++) seq_read_dealloc(&chroms.data[i]);
  readbuf_dealloc(&chroms);

  ctx_free(nw_scoring_flank);
  ctx_free(nw_scoring_allele);
  alignment_free(alignment);
  needleman_wunsch_free(nw_aligner);

  sam_close(samfh);
  free(bam_header);
  ctx_free(bam);

  strbuf_free(line);
  vcf_entry_dealloc(&invcf, num_samples);
  vcf_entry_dealloc(&outvcf, num_samples);

  return EXIT_SUCCESS;
}
