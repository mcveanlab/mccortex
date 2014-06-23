#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "seq_reader.h"
#include "call_file_reader.h"

#include "sam.h"
#include "seq-align/src/needleman_wunsch.h"

#define DEFAULT_MIN_MAPQ 30 /* min MAPQ considered */
#define DEFAULT_MAX_ALEN 500 /* max allele length */
#define DEFAULT_MAX_PDIFF 500 /* max allele length */

const char calls2vcf_usage[] =
"usage: "CMD" calls2vcf [options] <in.txt.gz> <ref.fa> [ref2.fa ...]\n"
"\n"
"  Convert a bubble or breakpoint call file to VCF. If input is a bubble file\n"
"  the --mapped <flanks.sam> argument is required.\n"
"\n"
"  -h, --help             This help message\n"
"  -o, --out <out.txt>    Save output graph file [default: STDOUT]\n"
"\n"
"  -f, --flanks <in.bam>  Mapped flanks in SAM or BAM file\n"
"  -Q, --min-mapq <Q>     Flank must map with MAPQ >= <Q> [default: "QUOTE_VALUE(DEFAULT_MIN_MAPQ)"]\n"
"  -A, --max-allele <M>   Max allele length considered [default: "QUOTE_VALUE(DEFAULT_MAX_ALEN)"]\n"
"  -D, --max-diff <D>     Max difference in path lengths [default: "QUOTE_VALUE(DEFAULT_MAX_PDIFF)"]\n"
"\n"
"  Alignment scoring:\n"
"  -m, --match <m>       [default:  1]\n"
"  -M, --mismatch <m>    [default: -2]\n"
"  -g, --gap-open <m>    [default: -4]\n"
"  -G, --gap-extend <m>  [default: -1]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
// command specific
  {"flanks",       required_argument, NULL, 'f'},
  {"min-mapq",     required_argument, NULL, 'Q'},
  {"max-allele",   required_argument, NULL, 'A'},
  {"max-diff",     required_argument, NULL, 'D'},
// alignment
  {"match",        required_argument, NULL, 'm'},
  {"mismatch",     required_argument, NULL, 'M'},
  {"gap-open",     required_argument, NULL, 'g'},
  {"gap-extend",   required_argument, NULL, 'G'},
  {NULL, 0, NULL, 0}
};

//
// Parameters
//
static const char *input_path = NULL;
static const char *out_path = NULL, default_out_path[] = "-";
// Filtering parameters
static size_t min_mapq = SIZE_MAX;
static size_t max_allele_len = SIZE_MAX;
static size_t max_path_diff = SIZE_MAX;
// Alignment parameters
static int nwmatch = 1, nwmismatch = -2, nwgapopen = -4, nwgapextend = -1;
// ref paths
static char **ref_paths;
static size_t num_ref_paths = 0;
// flank file
const char *sam_path = NULL;

//
// Reference genome
//
// Hash map of chromosome name -> sequence
static khash_t(ChromHash) *genome;
static ReadBuffer chroms;

// Input file format
static bool input_bubble_format = false; // false => breakpoint format

// Flank mapping
static samFile *samfh;
static bam_hdr_t *bam_header;
static bam1_t *bam;

// nw alignment
static nw_aligner_t *nw_aligner;
static alignment_t *alignment;
static scoring_t nw_scoring_flank, nw_scoring_allele;

// Temporary memory
// static StrBuf endflank;

//
// Statistics
//
// VCF printing stats
static size_t num_entries_read = 0, num_vars_printed = 0;
// Alignment stats
// static size_t num_nw_flank = 0, num_nw_allele = 0;

static void parse_cmdline_args(int argc, char **argv)
{
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
      case 'o':
        if(out_path != NULL) cmd_print_usage(NULL);
        out_path = optarg;
        break;
      case 'f': if(sam_path){die("%s set twice", cmd);} sam_path=optarg; break;
      case 'Q':
        if(min_mapq != SIZE_MAX) die("%s set twice", cmd);
        min_mapq = cmd_parse_arg_uint32(cmd, optarg);
        break;
      case 'A':
        if(max_allele_len != SIZE_MAX) die("%s set twice", cmd);
        max_allele_len = cmd_parse_arg_uint32(cmd, optarg);
        break;
      case 'D':
        if(max_path_diff != SIZE_MAX) die("%s set twice", cmd);
        max_path_diff = cmd_parse_arg_uint32(cmd, optarg);
        break;
      case 'g': nwmatch = cmd_parse_arg_int32(cmd, optarg); break;
      case 'G': nwmismatch = cmd_parse_arg_int32(cmd, optarg); break;
      case 'm': nwgapopen = cmd_parse_arg_int32(cmd, optarg); break;
      case 'M': nwgapextend = cmd_parse_arg_int32(cmd, optarg); break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" bubbles -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = default_out_path;
  if(min_mapq == SIZE_MAX) min_mapq = DEFAULT_MIN_MAPQ;
  if(max_allele_len == SIZE_MAX) max_allele_len = DEFAULT_MAX_ALEN;
  if(max_path_diff == SIZE_MAX) max_path_diff = DEFAULT_MAX_PDIFF;

  if(optind+2 > argc)
    cmd_print_usage("Require <in.txt.gz> and at least one reference");

  input_path = argv[optind];
  optind++;
  ref_paths = argv + optind;
  num_ref_paths = argc - optind;
}

// Setup pairwise aligner
static void nw_aligner_setup()
{
  nw_aligner = needleman_wunsch_new();
  alignment = alignment_create(1024);
  scoring_init(&nw_scoring_flank, nwmatch, nwmismatch, nwgapopen, nwgapextend,
               true, true, 0, 0, 0, 0);
  scoring_init(&nw_scoring_allele, nwmatch, nwmismatch, nwgapopen, nwgapextend,
               false, false, 0, 0, 0, 0);
}

// Clean up pairwise aligner
static void nw_aligner_destroy()
{
  alignment_free(alignment);
  needleman_wunsch_free(nw_aligner);
}

// Returns 1 on success
//         0 if not mapped
//        -1 on error
static int sam_fetch_coords(const CallFileEntry *centry,
                            const char **chr_name,
                            size_t *start, size_t *end,
                            bool *fw_strand)
{
  (void)centry;
  (void)chr_name;
  (void)start;
  (void)end;
  (void)fw_strand;

  if(sam_read1(samfh, bam_header, bam) < 0) {
    warn("We've run out of SAM entries!");
    return -1;
  }

  int unmapped = bam->core.flag & BAM_FUNMAP;
  int low_mapq = bam->core.qual < min_mapq;

  if(unmapped || low_mapq) return 0;

  // look up chromosome
  *chr_name = bam_header->target_name[bam->core.tid];

  // Check ref name exists
  khiter_t k = kh_get(ChromHash, genome, *chr_name);
  if(k == kh_end(genome)) {
    warn("Cannot find chrom [%s]", *chr_name);
    return -1;
  }

  // const read_t *chr = kh_value(genome, k);

  // DEV 4: is this left or right flank
  //        assume complete left flank for now

  // DEV 3: Find other flank position


  return 1;
}

// Returns 1 on success
//         0 if not mapped
//        -1 on error
static int brkpnt_fetch_coords(const CallFileEntry *centry,
                               const char **chr_name,
                               size_t *start, size_t *end,
                               bool *fw_strand)
{
  (void)centry;
  (void)chr_name;
  (void)start;
  (void)end;
  (void)fw_strand;

  // Read for 5p, 3p mapping
  // DEV 1:
  return 1;
}

static void parse_entries(gzFile gzin, FILE *fout)
{
  (void)fout;

  CallFileEntry centry;
  call_file_entry_alloc(&centry);
  const char *chrom_name = NULL;
  size_t start = 0, end = 0;
  bool fw_strand = true, mapped = false;

  for(; call_file_read(gzin, input_path, &centry); num_entries_read++)
  {
    // Read a corresponding SAM entry
    if(sam_path) {
      int ret = sam_fetch_coords(&centry, &chrom_name, &start, &end, &fw_strand);
      if(ret < 0) break;
      mapped = (ret > 0);
    }
    else {
      mapped = brkpnt_fetch_coords(&centry, &chrom_name, &start, &end, &fw_strand);
    }

    if(mapped)
    {
      // DEV 2a: Check coords

      // DEV 2b: align
      
    }
  }

  call_file_entry_dealloc(&centry);
}

static void flanks_sam_open()
{
  if(!futil_path_has_extension(sam_path, ".bam") &&
     !futil_path_has_extension(sam_path, ".sam"))
  {
    cmd_print_usage("Mapped flanks is not .sam or .bam file: %s", sam_path);
  }

  bool isbam = futil_path_has_extension(sam_path, ".bam");

  samfh = sam_open(sam_path, isbam ? "rb" : "rs");
  if(samfh == NULL) die("Cannot open SAM/BAM %s", sam_path);

  // Load BAM header
  bam_header = sam_hdr_read(samfh);
  bam = bam_init1();
}

static void flanks_sam_close()
{
  sam_close(samfh);
  free(bam_header);
  free(bam);
}

#define strncasecmp2(str1,str2) strncasecmp((str1),(str2),strlen(str2))
#define line_is_vcfhdr(str,tag) (strncasecmp2((str),(tag)) == 0)

static void read_input_header(gzFile gzin, StrBuf *hdr)
{
  call_file_read_hdr(gzin, input_path, hdr);

  // Check file format to ensure it is BubbleCTX or BreakpointCTX
  char *line0end = strchr(hdr->buff, '\n');
  if(line0end == NULL) die("Empty file header? [path: %s]", input_path);
  *line0end = '\0';

  if(!line_is_vcfhdr(hdr->buff,"##fileformat="))
    die("Expected ##fileformat= to be first line [path: %s]: %s", input_path, hdr->buff);

  if(strncasecmp2(hdr->buff+strlen("##fileformat="),"CtxBreakpoints") == 0)
    input_bubble_format = false;
  else if(strncasecmp2(hdr->buff+strlen("##fileformat="),"CtxBubbles") == 0)
    input_bubble_format = true;
  else
    die("Unknown input format [path: %s]: %s", input_path, hdr->buff);

  *line0end = '\n';
}

static void print_vcf_header(StrBuf *hdr, FILE *fout)
{
  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  fprintf(fout, "##fileformat=VCFv4.1\n##fileDate=%s\n", datestr);

  // Print header lines, stripping out fileformat= and fileDate=
  char *str = hdr->buff, *lineend = strchr(str, '\n');
  while(1)
  {
    if(lineend) *lineend = '\0';

    if(!line_is_vcfhdr(str,"##fileformat=") &&
       !line_is_vcfhdr(str,"##fileDate="))
    {
      fputs(str, fout);
      fputc('\n', fout);
    }

    if(!lineend) break;
    else { *lineend = '\n'; str = lineend+1; lineend = strchr(str, '\n'); }
  }

  fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n", fout);
}

// Check contig entries match reference
// input file (breakpoint / bubble file) header contains lines:
//   ##contig=<ID=blah,length=123>
// We check that these match the reference just loaded
static void brkpnt_check_refs_match(StrBuf *hdr)
{
  size_t reflen, chrlen;
  char *name, *len;
  size_t num_contigs_parsed = 0;

  // Print header lines, stripping out fileformat= and fileDate=
  char *str = hdr->buff, *lineend = strchr(str, '\n');
  while(1)
  {
    if(line_is_vcfhdr(str,"##contig=<ID="))
    {
      name = str + strlen("##contig=<ID=");
      len = strstr(name, ",length=");
      if(len) {
        *len = '\0';
        len += strlen(",length=");
        chrlen = strtol(len, NULL, 10);
        // Find chrom
        khiter_t k = kh_get(ChromHash, genome, name);
        if(k == kh_end(genome))
          die("Cannot find ref chrom: %s", name);
        else {
          reflen = kh_value(genome, k)->seq.end;
          if(reflen != chrlen) {
            die("Chrom lengths do not match %s input:%zu ref:%zu",
                name, chrlen, reflen);
          }
        }
      }
      else die("Cannot parse contig entry");
      num_contigs_parsed++;
    }

    if(!lineend) break;
    else { str = lineend+1; lineend = strchr(str, '\n'); }
  }

  if(num_contigs_parsed != chroms.len) {
    die("Number of chromosomes differ: %zu in header vs %zu in ref",
        num_contigs_parsed, chroms.len);
  }
}

int ctx_calls2vcf(int argc, char **argv)
{
  parse_cmdline_args(argc, argv);

  gzFile gzin = futil_gzopen_input(input_path);
  FILE *fout = futil_open_output(out_path);

  nw_aligner_setup();

  // Read file header
  StrBuf hdr;
  strbuf_alloc(&hdr, 2048);
  read_input_header(gzin, &hdr);

  status("Reading %s in %s format", futil_inpath_str(input_path),
         input_bubble_format ? "bubble" : "breakpoint");

  // Open flank file if it exists
  if(sam_path) flanks_sam_open();

  // Load reference genome
  readbuf_alloc(&chroms, 1024);
  genome = kh_init(ChromHash);
  seq_reader_load_ref_genome(ref_paths, num_ref_paths, &chroms, genome);

  if(!input_bubble_format) brkpnt_check_refs_match(&hdr);

  // Run
  print_vcf_header(&hdr, fout);
  parse_entries(gzin, fout);

  char num_vars_printed_str[50], num_entries_read_str[50];
  ulong_to_str(num_entries_read, num_entries_read_str);
  ulong_to_str(num_vars_printed, num_vars_printed_str);
  status("Read %s entries, printed %s vcf entries to: %s",
         num_entries_read_str, num_vars_printed_str, futil_outpath_str(out_path));

  // Finished - clean up
  gzclose(gzin);
  fclose(fout);
  readbuf_dealloc(&chroms);
  kh_destroy_ChromHash(genome);
  nw_aligner_destroy();

  if(sam_path) flanks_sam_close();

  // hide unused method warnings
  (void)kh_del_ChromHash;
  (void)kh_put_ChromHash;
  (void)kh_get_ChromHash;
  (void)kh_clear_ChromHash;
  (void)kh_destroy_ChromHash;
  (void)kh_init_ChromHash;

  return EXIT_SUCCESS;
}
