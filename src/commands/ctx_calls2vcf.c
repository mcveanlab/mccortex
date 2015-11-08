#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "seq_reader.h"
#include "call_file_reader.h"
#include "json_hdr.h"
#include "vcf_misc.h"

#include "aligned_call.h"
#include "decomp_breakpoint.h"
#include "decomp_bubble.h"

#include "htslib/sam.h"
#include "seq-align/src/needleman_wunsch.h"

#define DEFAULT_MIN_MAPQ 30 /* min MAPQ considered (bubble caller only) */
#define DEFAULT_MAX_ALIGN 500 /* max path/bubble_branch length */
#define DEFAULT_MAX_ALLELE 500 /* max ALT allele length */

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
"  -O, --out-fmt <f>    Format vcf|vcfgz|bcf|ubcf\n"
"\n"
"  -F, --flanks <in.bam>  Mapped flanks in SAM or BAM file (bubble caller only)\n"
"  -Q, --min-mapq <Q>     Flank must map with MAPQ >= <Q> [default: "QUOTE_VALUE(DEFAULT_MIN_MAPQ)"]\n"
"  -A, --max-align <M>    Max alignment attempted [default: "QUOTE_VALUE(DEFAULT_MAX_ALIGN)"]\n"
"  -L, --max-allele <M>   Max allele length printed [default: "QUOTE_VALUE(DEFAULT_MAX_ALLELE)"]\n"
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
  {"out-fmt",      required_argument, NULL, 'O'},
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"flanks",       required_argument, NULL, 'F'},
  {"min-mapq",     required_argument, NULL, 'Q'},
  {"max-align",    required_argument, NULL, 'A'},
  {"max-allele",   required_argument, NULL, 'L'},
  {"max-diff",     required_argument, NULL, 'D'},
// alignment
  {"match",        required_argument, NULL, 'm'},
  {"mismatch",     required_argument, NULL, 'M'},
  {"gap-open",     required_argument, NULL, 'g'},
  {"gap-extend",   required_argument, NULL, 'G'},
  {NULL, 0, NULL, 0}
};


static void print_acall_stats(const DecomposeStats *stats)
{
  char n0[50], n1[50];
  status("[aligned] Calls mapped:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->ncalls_mapped, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->ncalls_mapped, stats->ncalls));
  status("[aligned] Ref allele too long:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->ncalls_ref_allele_too_long, n0),
         ulong_to_str(stats->ncalls_mapped, n1),
         safe_percent(stats->ncalls_ref_allele_too_long, stats->ncalls_mapped));
  // Alternative alleles before being decomposed
  uint64_t ncalls_good = stats->ncalls_mapped - stats->ncalls_ref_allele_too_long;
  status("[aligned] Alt. alleles per call:  %s / %s (%6.2f)",
         ulong_to_str(stats->nlines, n0),
         ulong_to_str(ncalls_good, n1),
         safe_percent(stats->nlines, ncalls_good)/100);
  status("[aligned] Alt. alleles too long:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nlines_too_long, n0),
         ulong_to_str(stats->nlines, n1),
         safe_percent(stats->nlines_too_long, stats->nlines));
  status("[aligned] Alt. alleles match ref:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nlines_match_ref, n0),
         ulong_to_str(stats->nlines, n1),
         safe_percent(stats->nlines_match_ref, stats->nlines));
  status("[aligned] Alt. alleles mapped:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nlines_mapped, n0),
         ulong_to_str(stats->nlines, n1),
         safe_percent(stats->nlines_mapped, stats->nlines));
  // Decomposed variants
  status("[aligned] ALTs too long:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nallele_too_long, n0),
         ulong_to_str(stats->nvars, n1),
         safe_percent(stats->nallele_too_long, stats->nvars));
  status("[aligned] ALTs printed:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nvars_printed, n0),
         ulong_to_str(stats->nvars, n1),
         safe_percent(stats->nvars_printed, stats->nvars));
}

static void print_bubble_stats(const DecompBubbleStats *stats)
{
  char n0[50], n1[50];
  status("[bubbles] 5' unmapped:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflank5p_unmapped, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflank5p_unmapped, stats->ncalls));
  status("[bubbles] 5' low MAPQ:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflank5p_lowqual, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflank5p_lowqual, stats->ncalls));
  status("[bubbles] 3' multi-hits:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflank3p_multihits, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflank3p_multihits, stats->ncalls));
  status("[bubbles] 3' not found:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflank3p_not_found, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflank3p_not_found, stats->ncalls));
  status("[bubbles] flanks overlap too much:  %s / %s (%6.2f%%) (DEL from ref)",
         ulong_to_str(stats->nflanks_overlap_too_much, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflanks_overlap_too_much, stats->ncalls));
  // success
  status("[bubbles] bubbles mapped:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->ncalls_mapped, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->ncalls_mapped, stats->ncalls));
  status("[bubbles] 3' kmer exact match:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflank3p_exact_found, n0),
         ulong_to_str(stats->ncalls_mapped, n1),
         safe_percent(stats->nflank3p_exact_found, stats->ncalls_mapped));
  status("[bubbles] 3' alignment found:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflank3p_approx_found, n0),
         ulong_to_str(stats->ncalls_mapped, n1),
         safe_percent(stats->nflank3p_approx_found, stats->ncalls_mapped));
}

static void print_breakpoint_stats(const DecompBreakpointStats *stats)
{
  char n0[50], n1[50];
  status("[brkpnts] flanks not uniquely mapped:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflanks_not_uniquely_mapped, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflanks_not_uniquely_mapped, stats->ncalls));
  status("[brkpnts] flanks map to diff chroms:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflanks_diff_chroms, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflanks_diff_chroms, stats->ncalls));
  status("[brkpnts] flanks map to diff strands:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflanks_diff_strands, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflanks_diff_strands, stats->ncalls));
  status("[brkpnts] flanks overlap too much:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->nflanks_overlap_too_much, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->nflanks_overlap_too_much, stats->ncalls));
  // success
  status("[brkpnts] breakpoints mapped:  %s / %s (%6.2f%%)",
         ulong_to_str(stats->ncalls_mapped, n0),
         ulong_to_str(stats->ncalls, n1),
         safe_percent(stats->ncalls_mapped, stats->ncalls));
}


static bcf_hdr_t* make_vcf_hdr(cJSON *json, const char *in_path,
                               bool is_breakpoint, size_t kmer_size,
                               char const*const* ref_paths, size_t nref_paths,
                               read_t *chroms, size_t nchroms)
{
  ctx_assert(json != NULL);

  StrBuf hdrbuf;
  strbuf_alloc(&hdrbuf, 1024);

  char datestr[9];
  time_t date = time(NULL);
  strftime(datestr, 9, "%Y%m%d", localtime(&date));

  strbuf_append_str(&hdrbuf, "##fileformat=VCFv4.2\n##fileDate=");
  strbuf_append_str(&hdrbuf, datestr);
  strbuf_append_str(&hdrbuf, "\n");

  // Print commands used to generate header
  cJSON *commands = json_hdr_get(json, "commands", cJSON_Array, in_path);
  cJSON *command = commands->child;

  // Print this command
  char keystr[8];
  char *prevstr = NULL;
  size_t i;

  if(command) {
    cJSON *key = json_hdr_get(command, "key", cJSON_String, in_path);
    prevstr = key->valuestring;
  }

  // Print command entry for this command
  strbuf_append_str(&hdrbuf, "##mccortex_");
  strbuf_append_str(&hdrbuf, hex_rand_str(keystr, sizeof(keystr)));
  strbuf_append_str(&hdrbuf, "=<prev=\"");
  strbuf_append_str(&hdrbuf, prevstr ? prevstr : "NULL");
  strbuf_append_str(&hdrbuf, "\",cmd=\"");
  strbuf_append_str(&hdrbuf, cmd_get_cmdline());
  strbuf_append_str(&hdrbuf, "\",cwd=\"");
  strbuf_append_str(&hdrbuf, cmd_get_cwd());
  strbuf_append_str(&hdrbuf, "\",version="CTX_VERSION">\n");

  // Print previous commands
  vcf_hdrtxt_append_commands(command, &hdrbuf, in_path);

  // Print field definitions
  if(is_breakpoint)
    strbuf_append_str(&hdrbuf, "##INFO=<ID=BRKPNT,Number=1,Type=String,Description=\"Breakpoint call\">\n");
  else
    strbuf_append_str(&hdrbuf, "##INFO=<ID=BUBBLE,Number=1,Type=String,Description=\"Bubble call\">\n");

  strbuf_sprintf(&hdrbuf, "##INFO=<ID=K%zu,Number=0,Type=Flag,Description=\"Found at k=%zu\">\n", kmer_size, kmer_size);

  strbuf_append_str(&hdrbuf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  strbuf_append_str(&hdrbuf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");

  // Print reference paths
  strbuf_append_str(&hdrbuf, "##reference=");
  strbuf_append_str(&hdrbuf, ref_paths[0]);
  for(i = 1; i < nref_paths; i++) {
    strbuf_append_char(&hdrbuf, ',');
    strbuf_append_str(&hdrbuf, ref_paths[i]);
  }
  strbuf_append_str(&hdrbuf, "\n");

  // Print contigs lengths
  for(i = 0; i < nchroms; i++) {
    strbuf_sprintf(&hdrbuf, "##contig=<ID=%s,length=%zu>\n",
                   chroms[i].name.b, chroms[i].seq.end);
  }

  // Print VCF column header
  strbuf_append_str(&hdrbuf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

  if(is_breakpoint)
  {
    // Print a column for each sample
    cJSON *graph_json   = json_hdr_get(json,       "graph",   cJSON_Object, in_path);
    cJSON *colours_json = json_hdr_get(graph_json, "colours", cJSON_Array,  in_path);
    cJSON *colour_json  = colours_json->child;
    if(colour_json == NULL) die("Missing colours");
    for(; colour_json; colour_json = colour_json->next)
    {
      if(!json_hdr_colour_is_ref(colour_json)) {
        cJSON *sample_json = json_hdr_get(colour_json, "sample", cJSON_String, in_path);
        strbuf_append_str(&hdrbuf, "\t");
        strbuf_append_str(&hdrbuf, sample_json->valuestring);
      }
    }
  }

  strbuf_append_char(&hdrbuf, '\n');
  bcf_hdr_t *hdr = bcf_hdr_init("w");
  if(bcf_hdr_parse(hdr, hdrbuf.b) != 0) die("Cannot construct VCF header");

  strbuf_dealloc(&hdrbuf);

  return hdr;
}

// Check contig entries match reference
// We check that these match the reference just loaded
static void brkpnt_check_refs_match(cJSON *json,
                                    const ChromHash *genome,
                                    const char *path)
{
  cJSON *version = json_hdr_get(json, "format_version", cJSON_Number, path);
  if(version->valueint <= 2) return;

  cJSON *command = json_hdr_get_curr_cmd(json, path);
  cJSON *brkpnts = json_hdr_get(command, "breakpoints", cJSON_Object, path);
  cJSON *contigs = json_hdr_get(brkpnts, "contigs",     cJSON_Array,  path);
  cJSON *contig;
  size_t num_chroms = 0;

  for(contig = contigs->child; contig; contig = contig->next, num_chroms++)
  {
    cJSON *id  = json_hdr_get(contig, "id",     cJSON_String, path);
    cJSON *len = json_hdr_get(contig, "length", cJSON_Number, path);
    const read_t *r = chrom_hash_fetch(genome, id->valuestring);
    if(r->seq.end != (size_t)len->valueint) {
      die("Chrom lengths do not match %s input:%li ref:%zu",
          id->valuestring, len->valueint, r->seq.end);
    }
  }

  if(num_chroms != kh_size(genome)) {
    die("Number of chromosomes differ: %zu in header vs %zu in ref",
        num_chroms, (size_t)kh_size(genome));
  }
}

int ctx_calls2vcf(int argc, char **argv)
{
  const char *in_path = NULL, *out_path = NULL, *out_type = NULL;
  // Filtering parameters
  int32_t min_mapq = -1, max_align_len = -1, max_allele_len = -1;
  // Alignment parameters
  int nwmatch = 1, nwmismatch = -2, nwgapopen = -4, nwgapextend = -1;
  // ref paths
  char const*const* ref_paths = NULL;
  size_t nref_paths = 0;
  // flank file
  const char *sam_path = NULL;

  //
  // Things we figure out by looking at the input
  //
  bool isbubble = false;
  // samples in VCF, (0 for bubble, does not include ref in breakpoint calls)
  size_t i, kmer_size, num_samples;

  //
  // Reference genome
  //
  // Hash map of chromosome name -> sequence
  ChromHash *genome;
  ReadBuffer chroms;

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
      case 'O': cmd_check(!out_type, cmd); out_type = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'F': cmd_check(!sam_path,cmd); sam_path = optarg; break;
      case 'Q': cmd_check(min_mapq < 0,cmd); min_mapq = cmd_uint32(cmd, optarg); break;
      case 'A': cmd_check(max_align_len  < 0,cmd); max_align_len  = cmd_uint32(cmd, optarg); break;
      case 'L': cmd_check(max_allele_len < 0,cmd); max_allele_len = cmd_uint32(cmd, optarg); break;
      case 'm': nwmatch = cmd_int32(cmd, optarg); break;
      case 'M': nwmismatch = cmd_int32(cmd, optarg); break;
      case 'g': nwgapopen = cmd_int32(cmd, optarg); break;
      case 'G': nwgapextend = cmd_int32(cmd, optarg); break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" calls2vcf -h` for help. Bad option: %s", argv[optind-1]);
      default: ctx_assert2(0, "shouldn't reach here: %c", c);
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(max_align_len  < 0) max_align_len  = DEFAULT_MAX_ALIGN;
  if(max_allele_len < 0) max_allele_len = DEFAULT_MAX_ALLELE;

  if(optind+2 > argc)
    cmd_print_usage("Require <in.txt.gz> and at least one reference");

  in_path = argv[optind++];
  ref_paths = (char const*const*)argv + optind;
  nref_paths = argc - optind;

  // These functions call die() on error
  gzFile gzin = futil_gzopen(in_path, "r");

  // Read call file header
  cJSON *json = json_hdr_load(gzin, in_path);

  // Check we can handle the kmer size
  kmer_size = json_hdr_get_kmer_size(json, in_path);
  db_graph_check_kmer_size(kmer_size, in_path);

  // Get format (bubble or breakpoint file)
  cJSON *json_fmt = json_hdr_get(json, "file_format", cJSON_String, in_path);
  if(strcmp(json_fmt->valuestring,"CtxBreakpoints") == 0) isbubble = false;
  else if(strcmp(json_fmt->valuestring,"CtxBubbles") == 0) isbubble = true;
  else die("Unknown format: '%s'", json_fmt->valuestring);

  status("Reading %s in %s format", futil_inpath_str(in_path),
         isbubble ? "bubble" : "breakpoint");

  if(isbubble) {
    // bubble specific
    if(sam_path == NULL)
      cmd_print_usage("Require -F <flanks.sam> with bubble file");
    if(min_mapq < 0) min_mapq = DEFAULT_MIN_MAPQ;
  }
  else {
    // breakpoint specific
    if(min_mapq >= 0)
      cmd_print_usage("-Q,--min-mapq <Q> only valid with bubble calls");
  }

  // Open flank file if it exists
  htsFile *samfh = NULL;
  bam_hdr_t *bam_hdr = NULL;
  bam1_t *mflank = NULL;

  if(sam_path)
  {
    if((samfh = hts_open(sam_path, "r")) == NULL)
      die("Cannot open SAM/BAM %s", sam_path);

    // Load BAM header
    bam_hdr = sam_hdr_read(samfh);
    if(bam_hdr == NULL) die("Cannot load BAM header: %s", sam_path);
    mflank = bam_init1();
  }

  // Output VCF has 0 samples if bubbles file, otherwise has N where N is
  // number of samples/colours in the breakpoint graph
  size_t num_graph_samples = json_hdr_get_ncols(json, in_path);
  size_t num_graph_nonref = json_hdr_get_nonref_ncols(json, in_path);

  num_samples = 0;
  if(!isbubble) {
    // If last colour has "is_ref", drop number of samples by one
    num_samples = num_graph_nonref < num_graph_samples ? num_graph_samples-1
                                                       : num_graph_samples;
  }

  //
  // Open output file
  //
  if(!out_path) out_path = "-";
  int mode = vcf_misc_get_outtype(out_type, out_path);
  futil_create_output(out_path);
  htsFile *vcffh = hts_open(out_path, modes_htslib[mode]);

  status("[calls2vcf] Reading %s call file with %zu samples",
         isbubble ? "Bubble" : "Breakpoint", num_graph_samples);
  status("[calls2vcf] %zu sample output to: %s format: %s",
         num_samples, futil_outpath_str(out_path), hsmodes_htslib[mode]);

  // Load reference genome
  read_buf_alloc(&chroms, 1024);
  genome = chrom_hash_init();
  chrom_hash_load(ref_paths, nref_paths, &chroms, genome);

  // convert to upper case
  char *s;
  for(i = 0; i < chroms.len; i++)
    for(s = chroms.b[i].seq.b; *s; s++) *s = toupper(*s);

  if(!isbubble) brkpnt_check_refs_match(json, genome, in_path);

  bcf_hdr_t *vcfhdr = make_vcf_hdr(json, in_path, !isbubble, kmer_size,
                                   ref_paths, nref_paths,
                                   chroms.b, chroms.len);

  if(bcf_hdr_write(vcffh, vcfhdr) != 0) die("Cannot write VCF header");

  AlignedCall *call = acall_init();
  CallDecomp *aligner = call_decomp_init(vcffh, vcfhdr);

  scoring_t *scoring = call_decomp_get_scoring(aligner);
  scoring_init(scoring, nwmatch, nwmismatch, nwgapopen, nwgapextend,
               false, false, 0, 0, 0, 0);

  CallFileEntry centry;
  call_file_entry_alloc(&centry);

  char kmer_str[50];
  sprintf(kmer_str, ";K%zu", kmer_size);

  if(isbubble)
  {
    // Bubble calls
    DecompBubble *bubbles = decomp_bubble_init();

    // Set scoring for aligning 3' flank
    scoring = decomp_bubble_get_scoring(bubbles);
    scoring_init(scoring, nwmatch, nwmismatch, nwgapopen, nwgapextend,
                 true, true, 0, 0, 0, 0);

    while(call_file_read(gzin, in_path, &centry)) {
      do {
        if(sam_read1(samfh, bam_hdr, mflank) < 0)
          die("We've run out of SAM entries!");
      } while(mflank->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));

      // Align call
      strbuf_reset(&call->info);
      decomp_bubble_call(bubbles, genome, kmer_size, min_mapq,
                         &centry, mflank, bam_hdr, call);
      strbuf_append_str(&call->info, kmer_str);
      acall_decompose(aligner, call, max_align_len, max_allele_len);
    }

    // print bubble stats
    DecompBubbleStats *bub_stats = ctx_calloc(1, sizeof(*bub_stats));
    decomp_bubble_cpy_stats(bub_stats, bubbles);
    print_bubble_stats(bub_stats);
    ctx_free(bub_stats);

    decomp_bubble_destroy(bubbles);
  }
  else
  {
    // Breakpoint calls
    DecompBreakpoint *breakpoints = decomp_brkpt_init();

    while(call_file_read(gzin, in_path, &centry)) {
      strbuf_reset(&call->info);
      decomp_brkpt_call(breakpoints, genome, num_samples, &centry, call);
      strbuf_append_str(&call->info, kmer_str);
      acall_decompose(aligner, call, max_align_len, max_allele_len);
    }

    // print bubble stats
    DecompBreakpointStats *brk_stats = ctx_calloc(1, sizeof(*brk_stats));
    decomp_brkpt_cpy_stats(brk_stats, breakpoints);
    print_breakpoint_stats(brk_stats);
    ctx_free(brk_stats);

    decomp_brkpt_destroy(breakpoints);
  }

  // Print stats
  DecomposeStats *astats = ctx_calloc(1, sizeof(*astats));
  call_decomp_cpy_stats(astats, aligner);
  print_acall_stats(astats);
  ctx_free(astats);

  call_file_entry_dealloc(&centry);
  call_decomp_destroy(aligner);
  acall_destroy(call);

  // Finished - clean up
  cJSON_Delete(json);
  gzclose(gzin);

  bcf_hdr_destroy(vcfhdr);
  hts_close(vcffh);

  for(i = 0; i < chroms.len; i++) seq_read_dealloc(&chroms.b[i]);
  read_buf_dealloc(&chroms);
  chrom_hash_destroy(genome);

  if(sam_path) {
    hts_close(samfh);
    bam_hdr_destroy(bam_hdr);
    bam_destroy1(mflank);
  }

  return EXIT_SUCCESS;
}
