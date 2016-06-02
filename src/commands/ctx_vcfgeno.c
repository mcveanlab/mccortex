#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "vcf_misc.h"
#include "genotyping.h"

#include <math.h> // log()
#include <float.h> // DBL_MAX
#include "htslib/vcf.h"

// TODO:
// [x] cmdline: specify diploid / haploid chroms
// [x] cmdline: differentiate between kmer-coverage and bp-coverage
// [x] cmdline: specify sample coverage
// [x] cmdline: set seqn error rate
// [x] read VCF entries
// [x] initialise lookup tables for maths
// [x] apply models, output results
// [x] get tags, kmer size vcf header
// [x] use bcf_int32_vector_end to mix haploid / diploid
// [x] report GT likelihood
// [x] print basic statistics
// [x] print genotyping statistics
// [x] vcfcov: add kmer coverage and read length to VCF header
// [x] vcfgeno: get kmer coverage and mean read length from VCF header
// [ ] add population classifier
// [ ] report variant likelihood (QUAL column)

#define SUBCMD "vcfgeno"

#define DEFAULT_ERR_RATE "0.01"

const char vcfgeno_usage[] =
"usage: "CMD" "SUBCMD" [options] <in.vcf>\n"
"\n"
"  Genotype a VCF after running vcfcov to add sample coverage. VCF does not\n"
"  need to be sorted. Currently only supports biallelic and ploidy 1 or 2.\n"
"\n"
"  -h, --help           This help message\n"
"  -q, --quiet          Silence status output normally printed to STDERR\n"
"  -f, --force          Overwrite output files\n"
"  -o, --out <out.vcf>  Output file [default: STDOUT]\n"
"  -O, --out-fmt <f>    Format vcf|vcfgz|bcf|ubcf\n"
"  -E, --err <E>        List of sample error rates (per bp) ["DEFAULT_ERR_RATE"]\n"
"  -D, --cov <C>        List of genome coverage per colour\n"
"  -C, --kcov <C>       List of kmer coverage per colour\n"
"  -P, --ploidy <P>     <ploidy> or sample:chr:ploidy (can be used repeatedly)\n"
"                       'sample' and 'chr' can be comma-separated lists\n"
"                       '.' means all. Be careful: applied in order.\n"
"  -l, --llk            Print all log likelihoods\n"
"  -r, --rm-cov         Remove tags set by 'vcfcov' command\n"
"  -R, --read-len <R>   List to override read lengths [optional]\n"
"\n"
"  Notes: \n"
"  1. Lists are comma-separated. If given a single value it will apply to all colours.\n"
"  2. User must give exactly one of --cov, --kcov\n"
"  3. 'kmer coverage' is calculated by: D * (R - k + 1) / R  where:\n"
"       D is sequence depth (X);  R is read length;  k is kmer size\n"
"     Or calculate with: `"CMD" view <graph.ctx>`\n"
"     Or look at output from: "CMD" clean <graph.ctx>\n"
"\n"
" Human ploidy example:\n"
"   "CMD" "SUBCMD" \\\n"
"      --cov 21,42,30 \\\n"
"      --ploidy .:.:2 --ploidy .:Y:0 --ploidy John,Tom:X,Y:1 \\\n"
"      in.vcf > out.vcf\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"out-fmt",      required_argument, NULL, 'O'},
  {"force",        no_argument,       NULL, 'f'},
  {"err",          required_argument, NULL, 'E'},
  {"cov",          required_argument, NULL, 'D'},
  {"kcov",         required_argument, NULL, 'C'},
  {"ploidy",       required_argument, NULL, 'P'},
  {"read-len",     required_argument, NULL, 'R'},
  {"llk",          no_argument,       NULL, 'l'},
  {"rm-cov",       no_argument,       NULL, 'r'},
  {NULL, 0, NULL, 0}
};


// Tags to use
char kcovgs_ref_tag[10] = {0}, kcovgs_alt_tag[10] = {0};

// Statistics
uint64_t num_lines_read = 0, num_ALTs_read = 0;
uint64_t num_missing_tags = 0; // no GT fields
uint64_t num_missing_covgs = 0, num_non_biallelic = 0;
uint64_t ploidy_seen[3] = {0,0,0}, num_genotypes_printed = 0;


#define N_LFAC 1024
double *lnfac_table = NULL; // log natural of x factorial: ln(x!)

static void math_calcs_init()
{
  size_t i;
  // natural log of factorial
  lnfac_table = ctx_malloc(N_LFAC * sizeof(lnfac_table[0]));
  for(i = 0; i < N_LFAC; i++) lnfac_table[i] = lgamma(i+1);
}

static void math_calcs_destroy() {
  ctx_free(lnfac_table);
  lnfac_table = NULL;
}

#define lnfac(x) ((x) < N_LFAC ? lnfac_table[x] : lgamma((x)+1))

// logerr = log(seq_err)
static double llk_hom(uint64_t covg1, uint64_t covg2, double theta1, double logerr)
{
  double logtheta1 = log(theta1); // log(err*theta) = log(err)+log(theta)
  return covg1 * logtheta1 - theta1 - lnfac(covg1) + covg2 * (logerr + logtheta1);
}

static double llk_het(uint64_t covg1, uint64_t covg2, double theta1, double theta2)
{
  return covg1 * log(theta1/2) - theta1/2 - lnfac(covg1) +
         covg2 * log(theta2/2) - theta2/2 - lnfac(covg2);
}

// TODO: update for ploidy > 2, nalts > 2
// get number of possible genotypes given ploidy
// assumes ploidy <= 2
static inline size_t num_gts(size_t ploidy, size_t nalts)
{
  (void)nalts;
  return !ploidy ? 0 : (ploidy == 1 ? 2 : 3);
}

static inline void set_gts_missing(size_t ploidy, size_t nalts,
                                   int32_t *gts, size_t ngts,
                                   float *gllks, size_t ngllks,
                                   int32_t *gt_qual)
{
  size_t i;
  for(i = 0; i < ploidy; i++) gts[i] = bcf_gt_missing;
  if(ploidy < ngts) gts[ploidy] = bcf_int32_vector_end;
  *gt_qual = bcf_int32_missing;
  if(gllks) {
    size_t expgllks = num_gts(ploidy, nalts);
    for(i = 0; i < expgllks; i++) bcf_float_set_missing(gllks[i]);
    if(expgllks < ngllks) bcf_float_set_vector_end(gllks[expgllks]);
  }
}

/**
 * Genotype a single sample
 * param gts sample genotype goes into: gts[0..(ploidy-1)]
 * param rcovgs,acovgs are of length num_alts (fixed to 1 currently)
 * param ploidy is this samples ploidy on the current chromosome
 * param logerr is ln(sample_err_rate)
 * param readlenk is read length in kmers
 **/
static void genotype_sample_biallelic(bcf1_t *v, size_t sampleid,
                                      int32_t *gts, size_t ngts,
                                      float *gllks, size_t ngllks,
                                      int32_t *gt_qual,
                                      const int32_t *rcovgs,
                                      const int32_t *acovgs,
                                      double kcovg, uint8_t ploidy,
                                      size_t ksize, double logerr,
                                      size_t readlenk)
{
  ctx_assert(ploidy <= ngts);
  ctx_assert(v->n_allele == 2);

  if(rcovgs[0] == bcf_int32_missing || acovgs[0] == bcf_int32_missing)
  {
    num_missing_covgs++;
    set_gts_missing(ploidy, 2, gts, ngts, gllks, ngllks, gt_qual);
  }
  else if(ploidy == 0)
  {
    set_gts_missing(ploidy, 2, gts, ngts, gllks, ngllks, gt_qual);
  }
  else
  {
    if(!readlenk)
      die("Read length is zero for sample: %zu", sampleid);

    // Get rlen, alen in kmers
    size_t rlen = 0, alen = 0, rshift;
    rshift = trimmed_alt_lengths(v, 1, &rlen, &alen);

    uint64_t rlenk = hap_num_exp_kmers(v->pos+rshift, rlen, ksize);
    uint64_t alenk = hap_num_exp_kmers(v->pos+rshift, alen, ksize);

    // theta1 is expected number of reads arriving on ref allele
    // theta2 is expected number of reads arriving on alt allele
    double theta1, theta2;
    uint64_t rkcov, akcov;
    double llk[3]; // hom1, het, hom2
    int order[3] = {0,1,2};

    // convert kmer coverage to num. of reads arriving, est. read arrival rate
    theta1 = kcovg * rlenk / readlenk;
    theta2 = kcovg * alenk / readlenk;
    rkcov = rcovgs[0] * rlenk / readlenk;
    akcov = acovgs[0] * alenk / readlenk;

    // given: loc_b(x) = log_a(x) / log_a(b)
    // llk_*() return natural log, divide by ln(10)=M_LN10 to get in log10
    // log10 is required by the VCFv4.2 standard for the GL tag
    llk[0] = llk_hom(rkcov, akcov, theta1, logerr) / M_LN10;
    llk[1] = ploidy == 2 ? llk_het(rkcov, akcov, theta1, theta2) / M_LN10 : -DBL_MAX;
    llk[2] = llk_hom(akcov, rkcov, theta2, logerr) / M_LN10;

    if(llk[order[0]] > llk[order[1]]) SWAP(order[0], order[1]);
    if(llk[order[1]] > llk[order[2]]) SWAP(order[1], order[2]);
    if(llk[order[0]] > llk[order[1]]) SWAP(order[0], order[1]);

    // if haploid: g0, if diploid: g0/g1
    uint32_t g0 = (order[2] == 2), g1 = (order[2] > 0);

    // set GT quality to be difference between highest and second highest GT llk
    *gt_qual = (int32_t)(llk[order[2]] - llk[order[1]] + 0.5);

    gts[0] = bcf_gt_unphased(g0);
    if(ploidy == 2) gts[1] = bcf_gt_unphased(g1);
    if(ploidy < ngts) gts[ploidy] = bcf_int32_vector_end;

    // ndecplaces(x,100) gets to two decimal places
    if(gllks) {
      if(ploidy == 1) {
        gllks[0] = ndecplaces(llk[0],100);
        gllks[1] = ndecplaces(llk[2],100);
        if(ngllks > 2) bcf_float_set_vector_end(gllks[2]);
      } else {
        gllks[0] = ndecplaces(llk[0],100);
        gllks[1] = ndecplaces(llk[1],100);
        gllks[2] = ndecplaces(llk[2],100);
        if(ngllks > 3) bcf_float_set_vector_end(gllks[3]);
      }
    }

    num_genotypes_printed++; // non-missing genotype printed
  }
}

static void genotype_vcf(htsFile *vcffh, bcf_hdr_t *vcfhdr, htsFile *outfh,
                         const double *kcovgs, const double *log_errs,
                         const uint8_t *const* ploidy_mat,
                         size_t *readlensk,
                         size_t max_ploidy, size_t kmer_size,
                         bool add_gllks, bool rm_vcfcov_tags)
{
  size_t nalts, s, nsamples = bcf_hdr_nsamples(vcfhdr);
  bcf1_t *v = bcf_init();
  int a, b;

  // ref / alt kmer coverage, one per ALT allele
  int nkcovr = 0, nkcova = 0;
  int32_t *kcovr = NULL, *kcova = NULL;

  size_t max_gts = num_gts(max_ploidy, 2);
  size_t ngllks = nsamples * max_gts;

  int ngts = nsamples * max_ploidy;
  int32_t *gts = ctx_calloc(ngts, sizeof(gts[0]));
  int32_t *gtquals = ctx_calloc(nsamples, sizeof(gtquals[0]));
  float *gllks = add_gllks ? ctx_calloc(ngllks, sizeof(gllks[0])) : NULL;

  // Initialise lookup tables
  math_calcs_init();

  // read and print
  while(bcf_read(vcffh, vcfhdr, v) == 0)
  {
    bcf_unpack(v, BCF_UN_ALL);
    num_lines_read++;
    num_ALTs_read += v->n_allele-1;

    if(v->n_allele != 2) { num_non_biallelic++; continue; }
    nalts = v->n_allele-1;

    // TODO: if add_gllks and nalts > 1, may need to resize gllks

    a = bcf_get_format_int32(vcfhdr, v, kcovgs_ref_tag, &kcovr, &nkcovr);
    b = bcf_get_format_int32(vcfhdr, v, kcovgs_alt_tag, &kcova, &nkcova);
    if(a < 0 || b < 0) { num_missing_tags++; continue; }

    // loop over samples
    for(s = 0; s < nsamples; s++) {
      uint8_t ploidy = ploidy_mat[v->rid][s];
      ploidy_seen[ploidy]++;
      genotype_sample_biallelic(v, s,
                                gts+max_ploidy*s, max_ploidy,
                                gllks ? gllks+max_gts*s : NULL, max_gts,
                                gtquals+s,
                                kcovr+nalts*s, kcova+nalts*s,
                                kcovgs[s], ploidy, kmer_size, log_errs[s],
                                readlensk[s]);
    }

    // Update GTs and write out
    if(rm_vcfcov_tags) {
      // remove coverage tags added by vcfcov
      a = bcf_update_format_int32(vcfhdr, v, kcovgs_ref_tag, NULL, 0);
      b = bcf_update_format_int32(vcfhdr, v, kcovgs_alt_tag, NULL, 0);
    }
    if(bcf_update_genotypes(vcfhdr, v, gts, ngts) < 0)
      die("Cannot update GTs");
    if(bcf_update_format_int32(vcfhdr, v, "GQ", gtquals, nsamples) < 0)
      die("Cannot update GQs");
    if(gllks && bcf_update_format_float(vcfhdr, v, "GL", gllks, ngllks) < 0)
      die("Cannot update GLs");
    if(bcf_write(outfh, vcfhdr, v) != 0) die("Cannot write record");
  }

  free(kcovr);
  free(kcova);
  ctx_free(gts);
  ctx_free(gtquals);
  ctx_free(gllks);
  bcf_destroy(v);
}

static int match_list(const char *s, char const*const* list, size_t n)
{
  size_t i;
  for(i = 0; i < n; i++)
    if(strcmp(s,list[i]) == 0)
      return i;
  return -1;
}

// returns true on success, false if cannot parse
// Format <samples>:<chrs>:<ploidy> or <ploidy>
//   '<ploidy>' is equivalent to '.:.:<ploidy>'
// 'samples' and 'chrs' can be comma-separated lists or '.' (means ALL)
// Valid values:
//   1
//   2
//   John,Jane:chr1,chr2:2
//   John,Jane:.:2
//   .:.:2
//   John:X:1
//   John:Y:1
static bool parse_ploidy_arg(char *str, uint8_t **ploidy_mat,
                             char const*const* seqnames, size_t nseqs,
                             char const*const* samples, size_t nsamples)
{
  char *f[3], *end0 = NULL, *end1 = NULL, *sample, *chr;
  int s, s_beg, s_end, c, c_beg, c_end;
  size_t i, j, ploidy = 0;

  ctx_assert(str != NULL);

  // Single number (e.g. 2) on its own is equivalent to .:.:2
  if(parse_entire_size(str, &ploidy))
  {
    if(ploidy > 2) die("ploidy >2 not currently supported.");
    for(i = 0; i < nseqs; i++)
      for(j = 0; j < nsamples; j++)
        ploidy_mat[i][j] = ploidy;
    return true;
  }

  if((f[0] = strtok_r(str, ":", &end0)) == NULL) return false;
  if((f[1] = strtok_r(NULL,":", &end0)) == NULL) return false;
  if((f[2] = strtok_r(NULL,":", &end0)) == NULL) return false;

  // Parse ploidy
  if(!parse_entire_size(f[2], &ploidy)) return false;
  if(ploidy > 2) die("ploidy >2 not currently supported.");

  // Loop over chromosomes, then samples
  // both are comma separated arrays
  // '.' means ALL
  for(end0 = NULL; (chr = strtok_r(f[1], ",", &end0)); f[1] = NULL) {
    if(strcmp(chr,".") == 0) { c_beg = 0; c_end = nseqs; }
    else {
      c_beg = match_list(chr, seqnames, nseqs);
      if(c_beg < 0) die("Cannot find chrom: %s", chr);
      c_end = c_beg+1;
    }
    for(end1 = NULL; (sample = strtok_r(f[0], ",", &end1)); f[0] = NULL) {
      if(strcmp(sample,".") == 0) { s_beg = 0; s_end = nsamples; }
      else {
        s_beg = match_list(sample, samples, nsamples);
        if(s_beg < 0) die("Cannot find sample: %s", sample);
        s_end = s_beg+1;
      }
      for(c = c_beg; c < c_end; c++)
        for(s = s_beg; s < s_end; s++)
          ploidy_mat[c][s] = ploidy;
    }
  }
  return true;
}

// Given string of e.g. 'K31R' or 'K31A' (regex /^K(\d+)[RA]$/)
// return kmer size, otherwise return -1 for any other string
static int kmer_from_tag(const char *key)
{
  char *end;
  unsigned long kmer;
  if(key[0] != 'K') return -1;
  kmer = strtoul(key+1, &end, 10);
  if(!end || end == key+1 || (strcmp(end,"R") && strcmp(end,"A"))) return -1;
  return (kmer > (unsigned long)INT_MAX || (kmer&1) == 0 ? -1 : (int)kmer);
}

// find K31R, K31A tags in header
// Complain if there is more than one kmer present
static size_t kmer_from_hdr(const bcf_hdr_t *hdr)
{
  int i, j, n = 0, kmer, kmers[3];
  char *key, *val, *tags[3];

  // Read up to three tags matching: ^K\d+[AR]$ e.g. K31R
  for(i = 0; i < hdr->nhrec && n < 3; i++) {
    for(j = 0; j < hdr->hrec[i]->nkeys && n < 3; j++) {
      key = hdr->hrec[i]->keys[j];
      val = hdr->hrec[i]->vals[j];
      if(strcmp(key,"ID") == 0 && (kmer = kmer_from_tag(val)) >= 0) {
        kmers[n] = kmer;
        tags[n] = val;
        n++;
      }
    }
  }

  if(n == 0) die("VCF hdr tags missing (e.g. K31R, K31A for k=31): run vcfcov!");
  if(n > 2) die("Multiple kmer tags in hdr: %s,%s,%s",tags[0],tags[1],tags[2]);
  if(kmers[0] != kmers[1]) die("Conflicting tags: %s %s", tags[0], tags[1]);
  if(!strcmp(tags[0],tags[1])) die("duplicate tags?"); // htslib shouldn't allow

  return kmers[0];
}

// Returns true if mean read lengths found for all samples
// otherwise return false
bool readlens_from_hdr(const bcf_hdr_t *hdr, size_t *rls, const char *rltag)
{
  bcf_hrec_t *hrec;
  int keyidx;
  size_t i, nsamples = bcf_hdr_nsamples(hdr);
  bool missing = false;

  for(i = 0; i < nsamples; i++) {
    hrec = bcf_hdr_get_hrec(hdr, BCF_HL_STR, "ID", hdr->samples[i], "SAMPLE");
    if(!hrec) { missing = true; continue; }
    keyidx = bcf_hrec_find_key(hrec, rltag);
    if(keyidx < 0) { missing = true; continue; }
    if(!parse_entire_size(hrec->vals[keyidx], &rls[i])) {
      die("Sample '%s' has an invalid read length in header: %s",
          hdr->samples[i], hrec->vals[keyidx]);
    }
  }

  return !missing;
}

int ctx_vcfgeno(int argc, char **argv)
{
  const char *out_path = NULL, *out_type = NULL;
  const char *kcov_arg = NULL, *cov_arg = NULL;
  const char *err_arg = NULL, *readlen_arg = NULL;
  char *pl_args[argc];
  size_t npl_args = 0;
  bool add_gllks = false, rm_vcfcov_tags = false;

  // These are set after we have initially looped over args
  size_t *read_lens = NULL;
  double *err_rates = NULL; // sample error rate array
  double *kcovgs = NULL; // sample kmer coverage array
  uint8_t **ploidy_mat = NULL; // ploidy martix: [chrid][sample]

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;
  size_t i;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'O': cmd_check(!out_type, cmd); out_type = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'E': cmd_check(!err_arg, cmd); err_arg = optarg; break;
      case 'C': cmd_check(!kcov_arg, cmd); kcov_arg = optarg; break;
      case 'D': cmd_check(!cov_arg, cmd); cov_arg = optarg; break;
      case 'P': pl_args[npl_args++] = optarg; break;
      case 'R': cmd_check(!readlen_arg, cmd); readlen_arg = optarg; break;
      case 'l': cmd_check(!add_gllks, cmd); add_gllks = true; break;
      case 'r': cmd_check(!rm_vcfcov_tags, cmd); rm_vcfcov_tags = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" "SUBCMD" -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  char default_err[] = DEFAULT_ERR_RATE;

  // if(!err_arg) cmd_print_usage("Require '--err 0.01,0.005,...' argument");
  if(!err_arg) err_arg = default_err;
  if(!npl_args) cmd_print_usage("Require '--ploidy sample:chr:ploidy' argment");
  if(!kcov_arg && !cov_arg) cmd_print_usage("Require --kcov or --cov argument");
  if(optind+1 != argc) cmd_print_usage("Need to pass a single VCF");
  const char *inpath = argv[optind];

  size_t s, max_ploidy = 0, kmer_size = 0;

  // Open input VCF file
  htsFile *vcffh = hts_open(inpath, "r");
  if(vcffh == NULL) die("Cannot open VCF file: %s", inpath);
  bcf_hdr_t *vcfhdr = bcf_hdr_read(vcffh);
  if(vcfhdr == NULL) die("Cannot read VCF header: %s", inpath);
  size_t nsamples = bcf_hdr_nsamples(vcfhdr);

  if(nsamples == 0) die("No samples in VCF");
  status("[vcfgeno] VCF has %zu samples", nsamples);

  // Get reference contig names
  int r, n, nseqs = 0;
  const char **seqnames = bcf_hdr_seqnames(vcfhdr, &nseqs);

  // Parse read lengths if passed, otherwise read from header
  if(readlen_arg)
  {
    read_lens = ctx_calloc(nsamples, sizeof(read_lens[0]));
    n = parse_list_sizes(read_lens, nsamples, readlen_arg);
    if(n < 0 || (n != 1 && n != (int)nsamples)) {
      die("Bad --read-lens arg: needs 1 or %zu comma separated numbers: %s",
          nsamples, readlen_arg);
    }
    while(n < (int)nsamples) read_lens[n++] = read_lens[0];
  }

  // Parse error arg
  err_rates = ctx_calloc(nsamples, sizeof(err_rates[0]));
  n = parse_list_doubles(err_rates, nsamples, err_arg);
  if(n < 0 || (n != 1 && n != (int)nsamples)) {
    die("Bad --err arg: needs 1 or %zu comma separated numbers: %s",
        nsamples, err_arg);
  }
  while(n < (int)nsamples) err_rates[n++] = err_rates[0];

  // Calculate log(err_rate) for each sample
  double *log_errs = ctx_calloc(nsamples, sizeof(log_errs[0]));
  for(i = 0; i < nsamples; i++) log_errs[i] = log(err_rates[i]);

  // Initialise sample coverage array
  kcovgs = ctx_calloc(nsamples, sizeof(kcovgs[0]));
  n = parse_list_doubles(kcovgs, nsamples, cov_arg ? cov_arg : kcov_arg);
  if(n < 0 || (n != 1 && n != (int)nsamples)) {
    die("Bad %s arg: needs 1 or %zu comma separated numbers: %s",
        cov_arg ? "--cov" : "--kcov", nsamples, kcov_arg);
  }
  while(n < (int)nsamples) kcovgs[n++] = kcovgs[0];

  // Initialise ploidy martix
  ploidy_mat = ctx_calloc(nseqs, sizeof(ploidy_mat[0]));
  for(r = 0; r < nseqs; r++)
    ploidy_mat[r] = ctx_calloc(nsamples, sizeof(ploidy_mat[0][0]));

  StrBuf tmpstr = {0,0,0};
  strbuf_alloc(&tmpstr, 256);

  strbuf_reset(&tmpstr);
  for(i = 0; i < npl_args; i++) {
    strbuf_set(&tmpstr, pl_args[i]);
    if(!parse_ploidy_arg(pl_args[i], ploidy_mat,
                         (char const*const*)seqnames, nseqs,
                         (char const*const*)vcfhdr->samples, nsamples))
    {
      die("Poorly formatted '--ploidy %s'", tmpstr.b);
    }
  }

  // set max ploidy
  bool zero_ploidy_seen = false;
  for(r = 0; r < nseqs; r++) {
    uint8_t seq_ploidy_merge = true;
    for(s = 0; s < nsamples; s++) {
      max_ploidy = MAX2(max_ploidy, ploidy_mat[r][s]);
      seq_ploidy_merge |= ploidy_mat[r][s];
    }
    zero_ploidy_seen |= (seq_ploidy_merge == 0);
  }

  // Warn if some contigs have zero ploidy for everyone
  if(zero_ploidy_seen) warn("Some chromosomes have zero ploidy on all samples");
  if(max_ploidy == 0) die("ploidy is all zero");
  if(max_ploidy > 2) die("some ploidy is crazy high (>2)");

  // Get coverage tags from header and use to get kmer size
  kmer_size = kmer_from_hdr(vcfhdr);

  // make read length tag
  char readlen_tag[] = "mean_read_length";

  // Get read lengths from VCF header if not set
  size_t *hdr_readlens = ctx_calloc(nsamples, sizeof(hdr_readlens[0]));
  if(!readlens_from_hdr(vcfhdr, hdr_readlens, readlen_tag))
    status("[vcfgeno] Note: vcfcov SAMPLE= headers not found for all samples");
  if(!read_lens) read_lens = hdr_readlens;
  else {
    strbuf_reset(&tmpstr);
    for(i = 0; i < nsamples; i++) strbuf_sprintf(&tmpstr, "\t%zu", hdr_readlens[i]);
    status("[vcfgeno] Read lengths from VCF header:%s", tmpstr.b);
    ctx_free(hdr_readlens);
  }

  #define readklen(len,ks) (MAX2((len), (ks)) - (ks) + 1)

  // Convert genome coverage to kmer coverage,
  // or set to zero if we don't have read length
  if(cov_arg) {
    for(i = 0; i < nsamples; i++) {
      kcovgs[i] = !read_lens[i] ? 0 : kcovgs[i] * readklen(read_lens[i],kmer_size)
                                      / (double)read_lens[i];
    }
  }

  // convert read lengths to kmers
  // leave unset (0) read lengths as zero
  // otherwise round up to 1
  for(i = 0; i < nsamples; i++)
    if(read_lens[i])
      read_lens[i] = readklen(read_lens[i],kmer_size);

  // Reconstruct tags from kmersize
  sprintf(kcovgs_ref_tag, "K%zuR", kmer_size);
  sprintf(kcovgs_alt_tag, "K%zuA", kmer_size);

  vcf_misc_hdr_add_cmd(vcfhdr, cmd_get_cmdline(), cmd_get_cwd());

  bcf_hdr_append(vcfhdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(vcfhdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality: difference between highest and next most likely log10 genotype likelihoods\">");

  if(add_gllks)
    bcf_hdr_append(vcfhdr, "##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype log10 likelihoods for 0/0, 0/1, 1/1\">");

  if(rm_vcfcov_tags) {
    bcf_hdr_remove(vcfhdr, BCF_HL_FMT, kcovgs_ref_tag);
    bcf_hdr_remove(vcfhdr, BCF_HL_FMT, kcovgs_alt_tag);
  }

  //
  // Print details to user
  //
  // ploidy matrix
  if(ctx_msg_out) {
    pthread_mutex_lock(&ctx_biglock);
    timestamp();
    message("[vcfgeno] Ploidy matrix:\n");
    for(s = 0; s < nsamples; s++) message("\t%s", vcfhdr->samples[s]);
    message("\n");
    for(r = 0; r < nseqs; r++) {
      message("%s", seqnames[r]);
      for(s = 0; s < nsamples; s++) message("\t%u", ploidy_mat[r][s]);
      message("\n");
    }
    pthread_mutex_unlock(&ctx_biglock);
  }

  // Print coverage tags we are using
  status("[vcfgeno] Using tags: %s, %s; kmer: %zu",
         kcovgs_ref_tag, kcovgs_alt_tag, kmer_size);

  strbuf_reset(&tmpstr);
  for(i = 0; i < nsamples; i++) {
    strbuf_append_char(&tmpstr, '\t');
    strbuf_append_str(&tmpstr, vcfhdr->samples[i]);
  }
  status("[vcfgeno] Sample names:%s", tmpstr.b);

  strbuf_reset(&tmpstr);
  for(i = 0; i < nsamples; i++) strbuf_sprintf(&tmpstr, "\t%zu", read_lens[i]);
  status("[vcfgeno] Read lengths:%s", tmpstr.b);

  strbuf_reset(&tmpstr);
  for(i = 0; i < nsamples; i++) strbuf_sprintf(&tmpstr, "\t%.1f", kcovgs[i]);
  status("[vcfgeno] Kmer coverages:%s", tmpstr.b);

  strbuf_reset(&tmpstr);
  for(i = 0; i < nsamples; i++) strbuf_sprintf(&tmpstr, "\t%.4f", err_rates[i]);
  status("[vcfgeno] Seqn error rates: %s", tmpstr.b);

  status("[vcfgeno] max ploidy: %zu", max_ploidy);

  //
  // Open output file
  //
  if(!out_path) out_path = "-";
  int mode = vcf_misc_get_outtype(out_type, out_path);
  futil_create_output(out_path);
  htsFile *outfh = hts_open(out_path, modes_htslib[mode]);
  status("[vcfgeno] Output to: %s format: %s",
         futil_outpath_str(out_path), hsmodes_htslib[mode]);

  if(bcf_hdr_write(outfh, vcfhdr) != 0)
    die("Cannot write header to: %s", futil_outpath_str(out_path));

  // Ready to go
  genotype_vcf(vcffh, vcfhdr, outfh,
               kcovgs, log_errs, (const uint8_t *const*)ploidy_mat, read_lens,
               max_ploidy, kmer_size, add_gllks, rm_vcfcov_tags);

  uint64_t n_sample_alts = nsamples * num_ALTs_read;

  // Print statistics
  char n0[50], n1[50];
  status("[vcfgeno] Read %s VCF lines", ulong_to_str(num_lines_read, n0));
  status("[vcfgeno] Read %s ALT alleles", ulong_to_str(num_ALTs_read, n0));
  status("[vcfgeno] Lines skipped non-biallelic: %s / %s (%.2f%%)",
         ulong_to_str(num_non_biallelic, n0),
         ulong_to_str(num_lines_read, n1),
         safe_percent(num_non_biallelic, num_lines_read));
  status("[vcfgeno] Lines skipped no K..R/K..A fields: %s / %s (%.2f%%)",
         ulong_to_str(num_missing_tags, n0),
         ulong_to_str(num_lines_read, n1),
         safe_percent(num_missing_tags, num_lines_read));
  status("[vcfgeno] # missing kcov: %s / %s (%.2f%%)",
         ulong_to_str(num_missing_covgs, n0),
         ulong_to_str(n_sample_alts, n1),
         safe_percent(num_missing_covgs, n_sample_alts));
  status("[vcfgeno] # ploidy zero: %s / %s (%.2f%%)",
         ulong_to_str(ploidy_seen[0], n0),
         ulong_to_str(n_sample_alts, n1),
         safe_percent(ploidy_seen[0], n_sample_alts));
  status("[vcfgeno] # ploidy one: %s / %s (%.2f%%)",
         ulong_to_str(ploidy_seen[1], n0),
         ulong_to_str(n_sample_alts, n1),
         safe_percent(ploidy_seen[1], n_sample_alts));
  status("[vcfgeno] # ploidy two: %s / %s (%.2f%%)",
         ulong_to_str(ploidy_seen[2], n0),
         ulong_to_str(n_sample_alts, n1),
         safe_percent(ploidy_seen[2], n_sample_alts));
  status("[vcfgeno] # genotyped %s / %s (%.2f%%)",
         ulong_to_str(num_genotypes_printed, n0),
         ulong_to_str(n_sample_alts, n1),
         safe_percent(num_genotypes_printed, n_sample_alts));

  strbuf_dealloc(&tmpstr);
  ctx_free(read_lens);
  ctx_free(err_rates);
  ctx_free(log_errs);
  ctx_free(kcovgs);
  for(r = 0; r < nseqs; r++) ctx_free(ploidy_mat[r]);
  ctx_free(ploidy_mat);

  math_calcs_destroy();

  free(seqnames);
  bcf_hdr_destroy(vcfhdr);
  hts_close(vcffh);
  hts_close(outfh);

  return EXIT_SUCCESS;
}
