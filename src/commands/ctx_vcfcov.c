#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "gpath_checks.h"

#include "vcf_coverage.h"
#include "vcf_misc.h"

#include "htslib/vcf.h"
#include "htslib/faidx.h"

#define SUBCMD "vcfcov"

const char vcfcov_usage[] =
"usage: "CMD" "SUBCMD" [options] <in.vcf> <in.ctx> [in2.ctx ...]\n"
"\n"
"  Add coverage to a VCF using cortex graphs. It is recommended to use\n"
"  uncleaned graphs. The VCF must be sorted by position, with duplicates removed\n"
"  Indels ought to be left aligned to remove duplicates.\n"
"\n"
"  -h, --help             This help message\n"
"  -q, --quiet            Silence status output normally printed to STDERR\n"
"  -f, --force            Overwrite output files\n"
"  -m, --memory <mem>     Memory to use\n"
"  -n, --nkmers <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -o, --out <out.vcf>    Output file [default: STDOUT]\n"
"  -O, --out-fmt <f>      Format vcf|vcfgz|bcf|ubcf\n"
"  -r, --ref <ref.fa>     Reference file [required]\n"
"  -L, --max-var-len <A>  Only use alleles <= A bases long [default: "QUOTE_VALUE(DEFAULT_MAX_ALLELE_LEN)"]\n"
"  -N, --max-nvars <N>    Limit haplotypes to <= N variants [default: "QUOTE_VALUE(DEFAULT_MAX_GT_VARS)"]\n"
"  -M, --low-mem          Two-passes of VCF to only load needed kmers [default]\n"
"  -H, --high-mem         One-pass of VCF, all kmers loaded (when streaming VCF)\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"out-fmt",      required_argument, NULL, 'O'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"ref",          required_argument, NULL, 'r'},
  {"max-var-len",  required_argument, NULL, 'L'},
  {"max-nvars",    required_argument, NULL, 'N'},
  {"low-mem",      no_argument,       NULL, 'M'},
  {"high-mem",     no_argument,       NULL, 'H'},
  {NULL, 0, NULL, 0}
};

char kcov_ref_tag[10], kcov_alt_tag[10];

int ctx_vcfcov(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL, *out_type = NULL;

  uint32_t max_allele_len = 0, max_gt_vars = 0;
  char *ref_path = NULL;
  bool use_lowmem = false, use_himem = false;

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
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'r': cmd_check(!ref_path, cmd); ref_path = optarg; break;
      case 'L': cmd_check(!max_allele_len,cmd); max_allele_len = cmd_uint32(cmd,optarg); break;
      case 'N': cmd_check(!max_gt_vars,cmd); max_gt_vars = cmd_uint32(cmd,optarg); break;
      case 'M': cmd_check(!use_lowmem, cmd); use_lowmem = true; break;
      case 'H': cmd_check(!use_himem, cmd); use_himem = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" "SUBCMD" -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(ref_path == NULL) cmd_print_usage("Require a reference (-r,--ref <ref.fa>)");
  if(optind+2 > argc) cmd_print_usage("Require VCF and graph files");

  if(use_lowmem && use_himem)
    cmd_print_usage("Cannot use --low-mem and --high-mem together!");

  // Override number of kmers to use if --low-mem passed, since we calculate
  // number of kmers required anyway
  if(use_lowmem && memargs.num_kmers_set) {
    memargs.num_kmers_set = 0;
    memargs.num_kmers_set = false;
  }

  if(!max_allele_len) max_allele_len = DEFAULT_MAX_ALLELE_LEN;
  if(!max_gt_vars) max_gt_vars = DEFAULT_MAX_GT_VARS;

  status("[vcfcov] max allele length: %u; max number of variants: %u",
         max_allele_len, max_gt_vars);

  // open ref
  // index fasta with: samtools faidx ref.fa
  faidx_t *fai = fai_load(ref_path);
  if(fai == NULL) die("Cannot load ref index: %s / %s.fai", ref_path, ref_path);

  // Open input VCF file
  const char *vcf_path = argv[optind++];
  htsFile *vcffh = hts_open(vcf_path, "r");
  if(vcffh == NULL) die("Cannot open VCF file: %s", vcf_path);
  bcf_hdr_t *vcfhdr = bcf_hdr_read(vcffh);
  if(vcfhdr == NULL) die("Cannot read VCF header: %s", vcf_path);

  // default to low mem if we're not reading from STDIN
  bool low_mem = use_lowmem ||
                 (!use_himem && strcmp(vcf_path,"-"));

  // Test we can close and reopen files
  if(low_mem) {
    hts_close(vcffh);
    bcf_hdr_destroy(vcfhdr);
    if((vcffh = hts_open(vcf_path, "r")) == NULL)
      die("Cannot re-open VCF file: %s", vcf_path);
    if((vcfhdr = bcf_hdr_read(vcffh)) == NULL)
      die("Cannot re-read VCF header: %s", vcf_path);
  }

  //
  // Open graph files
  //
  const size_t num_gfiles = argc - optind;
  char **graph_paths = argv + optind;
  ctx_assert(num_gfiles > 0);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, NULL, 0, -1);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Covg)*8 * ncols;
  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        low_mem ? -1 : (int64_t)ctx_max_kmers,
                                        ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  // v=>vcf, z=>compressed vcf, b=>bcf, bu=>uncompressed bcf
  int mode = vcf_misc_get_outtype(out_type, out_path);
  futil_create_output(out_path);
  htsFile *outfh = hts_open(out_path, modes_htslib[mode]);
  status("[vcfcov] Output format: %s", hsmodes_htslib[mode]);


  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_COVGS);

  //
  // Set up tag names
  //

  // *R => ref, *A => alt
  sprintf(kcov_ref_tag, "K%zuR", db_graph.kmer_size); // mean coverage
  sprintf(kcov_alt_tag, "K%zuA", db_graph.kmer_size);

  // #SAMPLE=<ID=...,K29KCOV=...,K29NK=...,K29RLK>
  // - K29_kcov is empirical kmer coverage
  // - K29_nkmers is the number of kmers in the sample
  // - mean_read_length is the mean read length in bases
  char sample_kcov_tag[20], sample_nk_tag[20], sample_rlk_tag[20];
  sprintf(sample_kcov_tag, "K%zu_kcov", db_graph.kmer_size); // mean coverage
  sprintf(sample_nk_tag, "K%zu_nkmers", db_graph.kmer_size);
  sprintf(sample_rlk_tag, "mean_read_length");

  //
  // Load kmers if we are using --low-mem
  //

  VcfCovStats st;
  memset(&st, 0, sizeof(st));
  VcfCovPrefs prefs = {.kcov_ref_tag = kcov_ref_tag,
                       .kcov_alt_tag = kcov_alt_tag,
                       .max_allele_len = max_allele_len,
                       .max_gt_vars = max_gt_vars,
                       .load_kmers_only = false};

  if(low_mem)
  {
    status("[vcfcov] Loading kmers from VCF+ref");

    prefs.load_kmers_only = true;
    vcfcov_file(vcffh, vcfhdr, NULL, NULL, vcf_path, fai,
                NULL, &prefs, &st, &db_graph);

    // Close files
    hts_close(vcffh);
    bcf_hdr_destroy(vcfhdr);

    // Re-open files
    if((vcffh = hts_open(vcf_path, "r")) == NULL)
      die("Cannot re-open VCF file: %s", vcf_path);
    if((vcfhdr = bcf_hdr_read(vcffh)) == NULL)
      die("Cannot re-read VCF header: %s", vcf_path);

    prefs.load_kmers_only = false;
  }

  //
  // Load graphs
  //
  GraphLoadingStats gstats;
  memset(&gstats, 0, sizeof(gstats));

  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);
  gprefs.must_exist_in_graph = low_mem;

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &gstats);
    graph_file_close(&gfiles[i]);
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  //
  // Set up VCF header / graph matchup
  //
  size_t *samplehdrids = ctx_malloc(db_graph.num_of_cols * sizeof(size_t));

  // Add samples to vcf header
  bcf_hdr_t *outhdr = bcf_hdr_dup(vcfhdr);
  bcf_hrec_t *hrec;
  int sid;
  char hdrstr[200];

  for(i = 0; i < db_graph.num_of_cols; i++) {
    char *sname = db_graph.ginfo[i].sample_name.b;
    if((sid = bcf_hdr_id2int(outhdr, BCF_DT_SAMPLE, sname)) < 0) {
      bcf_hdr_add_sample(outhdr, sname);
      sid = bcf_hdr_id2int(outhdr, BCF_DT_SAMPLE, sname);
    }
    samplehdrids[i] = sid;

    // Add SAMPLE field
    hrec = bcf_hdr_get_hrec(outhdr, BCF_HL_STR, "ID", sname, "SAMPLE");

    if(hrec == NULL) {
      sprintf(hdrstr, "##SAMPLE=<ID=%s,%s=%"PRIu64",%s=%"PRIu64",%s=%zu>", sname,
              sample_kcov_tag,
              gstats.nkmers[i] ? gstats.sumcov[i] / gstats.nkmers[i] : 0,
              sample_nk_tag, gstats.nkmers[i],
              sample_rlk_tag, (size_t)db_graph.ginfo[i].mean_read_length);
      bcf_hdr_append(outhdr, hdrstr);
    }
    else {
      // mean kcovg
      sprintf(hdrstr, "%"PRIu64, gstats.sumcov[i] / gstats.nkmers[i]);
      vcf_misc_add_update_hrec(hrec, sample_kcov_tag, hdrstr);
      // num kmers
      sprintf(hdrstr, "%"PRIu64, gstats.nkmers[i]);
      vcf_misc_add_update_hrec(hrec, sample_nk_tag, hdrstr);
      // mean read length in kmers
      sprintf(hdrstr, "%zu", (size_t)db_graph.ginfo[i].mean_read_length);
      vcf_misc_add_update_hrec(hrec, sample_rlk_tag, hdrstr);
    }

    status("[vcfcov] Colour %zu: %s [VCF column %zu]", i, sname, samplehdrids[i]);
  }

  // Add genotype format fields
  // One field per alternative allele

  sprintf(hdrstr, "##FORMAT=<ID=%s,Number=A,Type=Integer,"
          "Description=\"Coverage on ref (k=%zu): sum(kmer_covs) / exp_num_kmers\">\n",
          kcov_ref_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, hdrstr);
  sprintf(hdrstr, "##FORMAT=<ID=%s,Number=A,Type=Integer,"
          "Description=\"Coverage on alt (k=%zu): sum(kmer_covs) / exp_num_kmers\">\n",
          kcov_alt_tag, db_graph.kmer_size);
  bcf_hdr_append(outhdr, hdrstr);

  bcf_hdr_set_version(outhdr, "VCFv4.2");

  // Add command string to header
  vcf_misc_hdr_add_cmd(outhdr, cmd_get_cmdline(), cmd_get_cwd());

  if(bcf_hdr_write(outfh, outhdr) != 0)
    die("Cannot write header to: %s", futil_outpath_str(out_path));

  status("[vcfcov] Reading %s and adding coverage", vcf_path);

  // Reset stats and get coverage
  memset(&st, 0, sizeof(st));

  vcfcov_file(vcffh, vcfhdr, outfh, outhdr, vcf_path, fai,
              samplehdrids, &prefs, &st, &db_graph);

  // Print statistics
  char ns0[50], ns1[50];
  status("[vcfcov] Read %s VCF lines", ulong_to_str(st.nvcf_lines, ns0));
  status("[vcfcov] Read %s ALTs", ulong_to_str(st.nalts_read, ns0));
  status("[vcfcov] Used %s kmers", ulong_to_str(st.ngt_kmers, ns0));
  status("[vcfcov] ALTs used: %s / %s (%.2f%%)",
         ulong_to_str(st.nalts_loaded, ns0), ulong_to_str(st.nalts_read, ns1),
         st.nalts_read ? (100.0*st.nalts_loaded) / st.nalts_read : 0.0);
  status("[vcfcov] ALTs too long (>%ubp): %s / %s (%.2f%%)", max_allele_len,
         ulong_to_str(st.nalts_too_long, ns0), ulong_to_str(st.nalts_read, ns1),
         st.nalts_read ? (100.0*st.nalts_too_long) / st.nalts_read : 0.0);
  status("[vcfcov] ALTs too dense (>%u within %zubp): %s / %s (%.2f%%)",
         max_gt_vars, db_graph.kmer_size,
         ulong_to_str(st.nalts_no_covg, ns0), ulong_to_str(st.nalts_read, ns1),
         st.nalts_read ? (100.0*st.nalts_no_covg) / st.nalts_read : 0.0);
  status("[vcfcov] ALTs printed with coverage: %s / %s (%.2f%%)",
         ulong_to_str(st.nalts_with_covg, ns0), ulong_to_str(st.nalts_read, ns1),
         st.nalts_read ? (100.0*st.nalts_with_covg) / st.nalts_read : 0.0);

  status("[vcfcov] Saved to: %s\n", out_path);

  ctx_free(samplehdrids);
  graph_loading_stats_destroy(&gstats);

  bcf_hdr_destroy(vcfhdr);
  bcf_hdr_destroy(outhdr);
  hts_close(vcffh);
  hts_close(outfh);
  fai_destroy(fai);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
