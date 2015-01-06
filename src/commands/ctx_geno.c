#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "seq_reader.h"
#include "gpath_checks.h"

const char geno_usage[] =
"usage: "CMD" geno [options] <in.vcf> <in.ctx> [in2.ctx ...]\n"
"\n"
"  Genotype a VCF using cortex graphs. VCF must be sorted by position. \n"
"  VCF must be a file, not piped in. It is recommended to use uncleaned graphs.\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -o, --out <bub.txt.gz>  Output file [default: STDOUT]\n"
"  -r, --ref <ref.fa>      Reference file\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"ref",          required_argument, NULL, 'r'},
  {NULL, 0, NULL, 0}
};

int ctx_geno(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;

  seq_file_t *tmp_seq_file;
  SeqFilePtrBuffer ref_buf;
  seq_file_ptr_buf_alloc(&ref_buf, 16);

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
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'r':
        if((tmp_seq_file = seq_open(optarg)) == NULL)
          die("Cannot read --seq file %s", optarg);
        seq_file_ptr_buf_add(&ref_buf, tmp_seq_file);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" bubbles -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";

  if(optind+2 >= argc) cmd_print_usage("Require VCF, ref and graph files");

  // Open VCF file
  const char *vcf_path = argv[optind++];
  gzFile vcf_file = futil_gzopen(vcf_path, "r");

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

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Covg) * ncols;
  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        -1, -1,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  FILE *fout = futil_open_create(out_path, "w");

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_COVGS);

  // Load references
  ReadBuffer chrom_buf;
  read_buf_alloc(&chrom_buf, 512);
  seq_load_all_reads(ref_buf.data, ref_buf.len, &chrom_buf);

  // Close reference files
  for(i = 0; i < ref_buf.len; i++) seq_close(ref_buf.data[i]);
  seq_file_ptr_buf_dealloc(&ref_buf);

  // TODO: Load kmers from VCF + ref
  

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = true,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Seek to the start of VCF file
  gzseek(vcf_file, SEEK_SET, 0);

  // genotype calls

  status("  saved to: %s\n", out_path);
  
  gzclose(vcf_file);
  fclose(fout);

  read_buf_dealloc(&chrom_buf);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
