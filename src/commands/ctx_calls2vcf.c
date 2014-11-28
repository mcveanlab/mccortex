#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "seq_reader.h"
#include "call_file_reader.h"
#include "json_hdr.h"
#include "chrom_pos_list.h" // Parse chromosome position lists

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
"  -q, --quiet            Silence status output normally printed to STDERR\n"
"  -f, --force            Overwrite output files\n"
"  -o, --out <out.txt>    Save output graph file [default: STDOUT]\n"
"\n"
"  -F, --flanks <in.bam>  Mapped flanks in SAM or BAM file\n"
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
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"flanks",       required_argument, NULL, 'F'},
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

bool input_bubble_format = false;

//
// Reference genome
//
// Hash map of chromosome name -> sequence
static khash_t(ChromHash) *genome;
static ReadBuffer chroms;

// Flank mapping
static samFile *samfh;
static bam_hdr_t *bam_header;
static bam1_t *bam;

// nw alignment
static nw_aligner_t *nw_aligner;
static alignment_t *aln;
static scoring_t nw_scoring_flank, nw_scoring_allele;

//
// Statistics
//
// VCF printing stats
static size_t num_entries_read = 0, num_vars_printed = 0;
// Alignment stats
// static size_t num_nw_flank = 0, num_nw_allele = 0;

static size_t num_flanks_not_uniquely_mapped = 0;
static size_t num_flanks_diff_chroms = 0;
static size_t num_flanks_diff_strands = 0;
static size_t num_flanks_overlap_too_large = 0;

static size_t num_entries_well_mapped = 0;

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
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'F': cmd_check(!sam_path,cmd); sam_path = optarg; break;
      case 'Q': cmd_check(min_mapq == SIZE_MAX,cmd); min_mapq = cmd_uint32(cmd, optarg); break;
      case 'A': cmd_check(max_allele_len == SIZE_MAX,cmd); max_allele_len = cmd_uint32(cmd, optarg); break;
      case 'D': cmd_check(max_path_diff == SIZE_MAX, cmd); max_path_diff = cmd_uint32(cmd, optarg); break;
      case 'm': nwmatch = cmd_int32(cmd, optarg); break;
      case 'M': nwmismatch = cmd_int32(cmd, optarg); break;
      case 'g': nwgapopen = cmd_int32(cmd, optarg); break;
      case 'G': nwgapextend = cmd_int32(cmd, optarg); break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" calls2vcf -h` for help. Bad option: %s", argv[optind-1]);
      default: ctx_assert2(0, "shouldn't reach here: %c", c);
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
  aln = alignment_create(1024);
  scoring_init(&nw_scoring_flank, nwmatch, nwmismatch, nwgapopen, nwgapextend,
               true, true, 0, 0, 0, 0);
  scoring_init(&nw_scoring_allele, nwmatch, nwmismatch, nwgapopen, nwgapextend,
               false, false, 0, 0, 0, 0);
}

// Clean up pairwise aligner
static void nw_aligner_destroy()
{
  alignment_free(aln);
  needleman_wunsch_free(nw_aligner);
}

static bool sam_fetch_coords(const CallFileEntry *centry,
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

  if(unmapped || low_mapq) return false;

  *chr_name = bam_header->target_name[bam->core.tid];

  // TODO is this left or right flank
  //      assume complete left flank for now

  // TODO Find other flank position


  return true;
}

/**
 * Fetch the largest match from a breakpoint call
 * @param is line to be parsed '>seqname ... chr=...'
 * @param buf is temporary buffer
 * @param flank is used to return result of largest match
 * @return 1 on success, 0 if not mapped. Calls die() on error
 */
static int brkpnt_fetch_first_match(const char *line, ChromPosBuffer *buf,
                                    ChromPosOffset *flank)
{
  char *list = strstr(line, " chr=");
  if(list == NULL) die("Cannot find flank position: %s", line);
  // Parse chr=seq0b:1-20:+:1,seq0a:2-20:+:2
  if(chrom_pos_list_parse(list+5, buf) < 0) die("Invalid positions: %s", line);
  return chrom_pos_list_get_largest(buf, flank);
}

static bool brkpnt_fetch_coords(const CallFileEntry *centry,
                                ChromPosBuffer *chrposbuf,
                                const char **chrom, size_t *start, size_t *end,
                                bool *fw_strand)
{
  ChromPosOffset flank5p, flank3p;
  size_t n;

  if((n = call_file_num_lines(centry)) < 6) die("Fewer than 6 lines: %zu", n);

  char *line0 = call_file_get_line(centry,0);
  char *line2 = call_file_get_line(centry,2);

  bool success = (brkpnt_fetch_first_match(line0, chrposbuf, &flank5p) &&
                  brkpnt_fetch_first_match(line2, chrposbuf, &flank3p));

  // Didn't map uniquely, with mismatching chromosomes or strands
  if(!success) { num_flanks_not_uniquely_mapped++; return false; }
  if(strcmp(flank5p.chrom,flank3p.chrom) != 0) { num_flanks_diff_chroms++; return false; }
  if(flank5p.fw_strand != flank3p.fw_strand) { num_flanks_diff_strands++; return false; }

  // Copy results. ChromPosOffset coords are 1-based.
  *chrom = flank5p.chrom;
  *fw_strand = flank5p.fw_strand;
  if(flank5p.fw_strand) { *start = flank5p.end+1; *end = flank3p.start; }
  else                  { *start = flank3p.end+1; *end = flank5p.start; }

  // Convert to 0-based coords
  (*start)--;
  (*end)--;

  return true;
}

static void strbuf_append_dna(StrBuf *buf, const char *src, size_t len,
                              bool fw_strand)
{
  strbuf_ensure_capacity(buf, buf->end + len);
  if(fw_strand) {
    strbuf_append_strn(buf, src, len);
  } else {
    dna_revcomp_str(buf->b+buf->end, src, len);
    buf->end += len;
    buf->b[buf->end] = '\0';
  }
}

static size_t align_get_start(const char *ref, const char *alt, size_t len,
                              size_t offset)
{
  size_t i;
  for(i = offset; i < len; i++) {
    if(ref[i] != alt[i]) return offset;
  }
  return len;
}

static size_t align_get_end(const char *ref, const char *alt, size_t len,
                            size_t offset)
{
  size_t i;
  for(i = offset; i < len; i++) {
    if(ref[i] == alt[i]) return offset;
  }
  return len;
}

static size_t align_get_len(const char *allele, size_t len)
{
  size_t i, nbases = 0;
  for(i = 0; i < len; i++)
    if(allele[i] != '-')
      nbases++;
  return nbases;
}

// Print allele with previous base and no deletions
// 'A--CG-T' with prev_base 'C' => print 'CACGT'
static void print_vcf_allele(int prev_base, const char *allele, size_t len,
                             FILE *fout)
{
  size_t i;
  if(prev_base > 0) fputc((char)prev_base, fout);
  for(i = 0; i < len; i++)
    if(allele[i] != '-')
      fputc(allele[i], fout);
}

// @param vcf_pos is 1-based
// @param prev_base is -1 if SNP otherwise previous base
static void print_vcf_entry(const char *chrom_name, size_t vcf_pos, int prev_base,
                            const char *ref, const char *alt,
                            size_t aligned_len, FILE *fout)
{
  // CHROM POS ID REF ALT QUAL FILTER INFO
  fprintf(fout, "%s\t%zu\tvar%zu\t", chrom_name, vcf_pos, num_vars_printed);
  print_vcf_allele(prev_base, ref, aligned_len, fout);
  fputc('\t', fout);
  print_vcf_allele(prev_base, alt, aligned_len, fout);
  fputs("\t.\tPASS\t", fout);
  fputs(input_bubble_format ? "BUBBLE" : "BRKPNT", fout);
  fputs("\tGT\n", fout);
  num_vars_printed++;
}

// @param ref_pos is 0-based here
static void align_biallelic(const char *ref, const char *alt, size_t aligned_len,
                            const read_t *chr, size_t ref_pos, FILE *fout)
{
  size_t i, start, end = 0;
  size_t ref_allele_len, alt_allele_len;
  int prev_base, vcf_pos;
  bool is_snp;

  while((start = align_get_start(ref, alt, aligned_len, end)) < aligned_len)
  {
    // Update ref offset
    for(i = end; i < start; i++)
      if(ref[i] != '-') ref_pos++;

    end = align_get_end(ref, alt, aligned_len, start);

    ref_allele_len = align_get_len(ref+start, end-start);
    alt_allele_len = align_get_len(alt+start, end-start);
    is_snp = ref_allele_len == 1 && alt_allele_len == 1;

    vcf_pos = ref_pos+1; // Convert to 1-based

    if(!is_snp) {
      prev_base = ref_pos > 0 ? chr->seq.b[ref_pos-1] : 'N';
      vcf_pos--;
    } else {
      prev_base = -1;
    }

    print_vcf_entry(chr->name.b, vcf_pos, prev_base,
                    ref+start, alt+start, end-start, fout);

    ref_pos += ref_allele_len;
  }
}

static void align_entry(CallFileEntry *centry,
                        const char *chrom_name, size_t ref_start, size_t ref_end,
                        bool fw_strand,
                        StrBuf *tmpbuf, FILE *fout)
{
  size_t i, ncpy = 0;
  size_t nlines = call_file_num_lines(centry);
  ctx_assert2(!(nlines&1) && nlines >= 6, "Too few lines: %zu", nlines);
  char *flank5p = call_file_get_line(centry,1);
  char *flank3p = call_file_get_line(centry,3);
  size_t flank5p_len = call_file_line_len(centry,1);
  size_t flank3p_len = call_file_line_len(centry,3);
  bool cpy_flnk_5p = false;

  if(ref_start > ref_end) {
    ncpy = ref_start - ref_end;
    if(ncpy > flank5p_len && ncpy > flank3p_len) {
      num_flanks_overlap_too_large++;
      // printf("Copy too much %zu > %zu %zu\n", ncpy, flank5p_len, flank3p_len);
      return; // can't align
    }
    cpy_flnk_5p = (ncpy > flank5p_len);
    if(fw_strand == cpy_flnk_5p) ref_start -= ncpy;
    else                         ref_end   += ncpy;
  }

  num_entries_well_mapped++;

  // Fetch chromosome
  khiter_t k = kh_get(ChromHash, genome, chrom_name);
  if(k == kh_end(genome)) die("Cannot find chrom [%s]", chrom_name);
  const read_t *chr = kh_value(genome, k);

  // If not fw strand, we need to flip each allele

  // Deal with alleles one at a time vs ref
  // First allele stored in line 5:
  //   >flank5p\n<seq>\n>flank3p\n<seq>\n>allele1\n<seq>
  char *line, *seq;
  size_t linelen, seqlen;

  for(i = 5; i < nlines; i+=2)
  {
    line = call_file_get_line(centry,i);
    linelen = call_file_line_len(centry,i);

    if(ncpy == 0 && fw_strand) {
      seq = line;
      seqlen = linelen;
    } else {
      strbuf_reset(tmpbuf);

      // printf(" ref_start: %zu ref_end: %zu\n", ref_start, ref_end);
      // printf(" ncpy: %zu fw_strand: %i cpy_flnk_5p: %i\n", ncpy, fw_strand, cpy_flnk_5p);

      if(ncpy > 0) {
        if(fw_strand && cpy_flnk_5p)
          strbuf_append_dna(tmpbuf, flank5p+flank5p_len-ncpy, ncpy, fw_strand);
        else if(!fw_strand && !cpy_flnk_5p)
          strbuf_append_dna(tmpbuf, flank3p, ncpy, fw_strand);
      }

      // Copy allele
      strbuf_append_dna(tmpbuf, line, linelen, fw_strand);

      if(ncpy > 0) {
        if(!fw_strand && cpy_flnk_5p)
          strbuf_append_dna(tmpbuf, flank5p+flank5p_len-ncpy, ncpy, fw_strand);
        else if(fw_strand && !cpy_flnk_5p)
          strbuf_append_dna(tmpbuf, flank3p, ncpy, fw_strand);
      }

      seq = tmpbuf->b;
      seqlen = tmpbuf->end;
    }

    // printf("%.*s vs %.*s\n", (int)(ref_end-ref_start), chr->seq.b + ref_start,
    //                          (int)seqlen, seq);

    // Align chrom and seq
    needleman_wunsch_align2(chr->seq.b + ref_start, seq,
                            ref_end-ref_start, seqlen,
                            &nw_scoring_allele, nw_aligner, aln);

    // Break into variants and print VCF
    align_biallelic(aln->result_a, aln->result_b, aln->length,
                    chr, ref_start, fout);
  }
}

static void parse_entries(gzFile gzin, FILE *fout)
{
  (void)fout;

  CallFileEntry centry;
  call_file_entry_alloc(&centry);

  ChromPosBuffer chrposbuf;
  chrompos_buf_alloc(&chrposbuf, 32);

  StrBuf tmpbuf;
  strbuf_alloc(&tmpbuf, 1024);

  const char *chrom_name;
  size_t ref_start, ref_end;
  bool mapped, fw_strand;

  for(; call_file_read(gzin, input_path, &centry); num_entries_read++)
  {
    // Read a corresponding SAM entry
    if(sam_path) {
      mapped = sam_fetch_coords(&centry,
                                &chrom_name, &ref_start, &ref_end, &fw_strand);
    }
    else {
      mapped = brkpnt_fetch_coords(&centry, &chrposbuf,
                                   &chrom_name, &ref_start, &ref_end, &fw_strand);
    }

    if(mapped)
      align_entry(&centry, chrom_name, ref_start, ref_end, fw_strand, &tmpbuf, fout);
  }

  call_file_entry_dealloc(&centry);
  chrompos_buf_dealloc(&chrposbuf);
  strbuf_dealloc(&tmpbuf);
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

static cJSON* read_input_header(gzFile gzin)
{
  cJSON *json;
  StrBuf hdrstr;
  strbuf_alloc(&hdrstr, 1024);
  json_hdr_read(NULL, gzin, input_path, &hdrstr);
  json = cJSON_Parse(hdrstr.b);
  if(json == NULL) die("Invalid JSON header: %s", input_path);

  // Check we can handle the kmer size
  size_t kmer_size = json_hdr_get_kmer_size(json, input_path);
  db_graph_check_kmer_size(kmer_size, input_path);

  strbuf_dealloc(&hdrstr);

  return json;
}

static void print_vcf_header(cJSON *json, FILE *fout)
{
  ctx_assert(json != NULL);

  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  fprintf(fout, "##fileFormat=VCFv4.1\n##fileDate=%s\n", datestr);

  // Print commands used to generate header
  cJSON *commands = cJSON_GetObjectItem(json, "commands");
  if(commands == NULL || commands->type != cJSON_Array)
    die("Missing 'commands' field in JSON header");

  cJSON *command = commands->child;

  // Print this command
  char keystr[8];
  char *prevstr = NULL;
  size_t i;

  if(command) {
    cJSON *key = cJSON_GetObjectItem(command, "key");
    if(key == NULL || key->type != cJSON_String) die("Invalid 'key' field");
    prevstr = key->valuestring;
  }

  // Print command entry for this command
  fprintf(fout, "##CMD=<key=\"%s\",prev=\"%s\",cmd=\"%s\",cwd=\"%s\">\n",
          hex_rand_str(keystr, sizeof(keystr)),
          prevstr ? prevstr : "NULL",
          cmd_get_cmdline(), cmd_get_cwd());

  // Print previous commands
  for(; command != NULL; command = command->next) {
    cJSON *key = json_hdr_get(command, "key", cJSON_String, input_path);
    cJSON *cmd = json_hdr_get(command, "cmd", cJSON_Array, input_path);
    cJSON *cwd = json_hdr_get(command, "cwd", cJSON_String, input_path);
    cJSON *prev = json_hdr_get(command, "prev", cJSON_Array, input_path);
    prev = prev->child; // result could be NULL
    if(prev && prev->type != cJSON_String) die("Invalid 'prev' field");
    fprintf(fout, "##CMD=<key=\"%s\",prev=\"%s", key->valuestring, prev ? prev->valuestring : "NULL");
    if(prev) {
      while((prev = prev->next) != NULL) fprintf(fout, ";%s", prev->valuestring);
    }
    fprintf(fout, "\",cmd=\"");
    for(i = 0, cmd = cmd->child; cmd; cmd = cmd->next, i++) {
      if(i > 0) fputc(' ', fout);
      fputs(cmd->valuestring, fout);
    }
    fprintf(fout, "\",cwd=\"%s\">\n", cwd->valuestring);
  }

  // Print contigs lengths
  for(i = 0; i < chroms.len; i++) {
    fprintf(fout, "##contig=<id=%s,length=%zu>\n",
            chroms.data[i].name.b, chroms.data[i].seq.end);
  }

  // Print VCF column header
  fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n", fout);
}

// Check contig entries match reference
// We check that these match the reference just loaded
static void brkpnt_check_refs_match(cJSON *json)
{
  cJSON *brkpnts = json_hdr_get(json,    "breakpoints", cJSON_Object, input_path);
  cJSON *contigs = json_hdr_get(brkpnts, "contigs",     cJSON_Array,  input_path);
  cJSON *contig;
  size_t num_chroms = 0;

  for(contig = contigs->child; contig; contig = contig->next, num_chroms++)
  {
    cJSON *id  = json_hdr_get(contig, "id",     cJSON_String, input_path);
    cJSON *len = json_hdr_get(contig, "length", cJSON_Number, input_path);

    const char *chrom_name = id->valuestring;
    long chrom_len = len->valueint;
    size_t reflen;

    khiter_t k = kh_get(ChromHash, genome, chrom_name);
    if(k == kh_end(genome))
      die("Cannot find ref chrom: %s", chrom_name);
    else {
      reflen = kh_value(genome, k)->seq.end;
      if(reflen != (size_t)chrom_len) {
        die("Chrom lengths do not match %s input:%li ref:%zu",
            chrom_name, chrom_len, reflen);
      }
    }
  }

  if(num_chroms != chroms.len) {
    die("Number of chromosomes differ: %zu in header vs %zu in ref",
        num_chroms, chroms.len);
  }
}

int ctx_calls2vcf(int argc, char **argv)
{
  parse_cmdline_args(argc, argv);

  // These functions call die() on error
  gzFile gzin = futil_gzopen(input_path, "r");

  nw_aligner_setup();

  // Read file header
  cJSON *json = read_input_header(gzin);

  // Get format (bubble or breakpoint file)
  cJSON *json_fmt = json_hdr_get(json, "file_format", cJSON_String, input_path);
  if(strcmp(json_fmt->valuestring,"CtxBreakpoints") == 0) input_bubble_format = false;
  else if(strcmp(json_fmt->valuestring,"CtxBubbles") == 0) input_bubble_format = true;
  else die("Unknown format: '%s'", json_fmt->valuestring);

  status("Reading %s in %s format", futil_inpath_str(input_path),
         input_bubble_format ? "bubble" : "breakpoint");

  if(input_bubble_format && sam_path == NULL)
    cmd_print_usage("Require -f <flanks.sam> with bubble file");

  // Open flank file if it exists
  if(sam_path) flanks_sam_open();

  // Open output file
  FILE *fout = futil_open_create(out_path, "w");

  // Load reference genome
  read_buf_alloc(&chroms, 1024);
  genome = kh_init(ChromHash);
  seq_reader_load_ref_genome(ref_paths, num_ref_paths, &chroms, genome);

  if(!input_bubble_format) brkpnt_check_refs_match(json);

  // Run
  print_vcf_header(json, fout);
  parse_entries(gzin, fout);

  char num_vars_printed_str[50], num_entries_read_str[50];
  ulong_to_str(num_entries_read, num_entries_read_str);
  ulong_to_str(num_vars_printed, num_vars_printed_str);
  status("Read %s entries, printed %s vcf entries to: %s",
         num_entries_read_str, num_vars_printed_str, futil_outpath_str(out_path));

  char num_nt_uniq_map_str[50], num_diff_chr_str[50], num_diff_strands_str[50];
  char num_fl_overlaps_big_str[50], num_well_mapped_str[50];

  ulong_to_str(num_flanks_not_uniquely_mapped, num_nt_uniq_map_str);
  ulong_to_str(num_flanks_diff_chroms, num_diff_chr_str);
  ulong_to_str(num_flanks_diff_strands, num_diff_strands_str);
  ulong_to_str(num_flanks_overlap_too_large, num_fl_overlaps_big_str);
  ulong_to_str(num_entries_well_mapped, num_well_mapped_str);

  status("   %s / %s (%.2f%%) flank pairs contain one flank not mapped uniquely",
         num_nt_uniq_map_str, num_entries_read_str,
         (100.0 * num_flanks_not_uniquely_mapped) / num_entries_read);
  status("   %s / %s (%.2f%%) flank pairs map to diff chroms",
         num_diff_chr_str, num_entries_read_str,
         (100.0 * num_flanks_diff_chroms) / num_entries_read);
  status("   %s / %s (%.2f%%) flank pairs map to diff strands",
         num_diff_strands_str, num_entries_read_str,
         (100.0 * num_flanks_diff_strands) / num_entries_read);
  status("   %s / %s (%.2f%%) flank pairs overlap too much",
         num_fl_overlaps_big_str, num_entries_read_str,
         (100.0 * num_flanks_overlap_too_large) / num_entries_read);
  status("   %s / %s (%.2f%%) flank pairs map well",
         num_well_mapped_str, num_entries_read_str,
         (100.0 * num_entries_well_mapped) / num_entries_read);

  // Finished - clean up
  cJSON_Delete(json);
  gzclose(gzin);
  fclose(fout);
  read_buf_dealloc(&chroms);
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
