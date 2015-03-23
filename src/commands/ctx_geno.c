#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "seq_reader.h"
#include "gpath_checks.h"
#include "seq_reader.h"
#include "genotyping.h"

#include "htslib/vcf.h"

// DEV: specify sample/chrom ploidy

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

typedef struct {
  // _kmers is the number of kmers unique to the allele
  // _sumcovg is the sum of coverages on those kmers
  // _medcovg is the median coverage of those kmers
  size_t refkmers, refsumcovg;// refmedcovg;
  size_t altkmers, altsumcovg;// altmedcovg;
} GenoCovg;

void bkey_get_covg(BinaryKmer bkey, uint64_t altref_bits,
                   size_t ntgts, GenoCovg *gts,
                   int colour, const dBGraph *db_graph)
{
  size_t i;
  dBNode node = db_graph_find(db_graph, bkey);

  if(node.key != HASH_NOT_FOUND) {
    Covg covg = colour >= 0 ? db_node_get_covg(db_graph, node.key, colour)
                            : db_node_sum_covg(db_graph, node.key);

    for(i = 0; i < ntgts; i++, altref_bits >>= 2) {
      if((altref_bits & 3) == 1) { gts[i].refkmers++; gts[i].refsumcovg += covg; }
      if((altref_bits & 3) == 2) { gts[i].altkmers++; gts[i].altsumcovg += covg; }
    }
  }
}

  // if(fwrite(bkey.b, 1, sizeof(bkey.b), fout) != sizeof(bkey.b))
  //   die("Cannot write");


void add_new_var(GenoVar *var, int32_t pos, const char *ref, const char *alt)
{
  while(*ref && *ref == *alt) { ref++; alt++; pos++; }
  var->pos = pos;
  var->alt = alt;
  var->altlen = strlen(alt);
  var->reflen = strlen(ref);
}

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

  if(optind+2 > argc) cmd_print_usage("Require VCF, ref and graph files");

  // Open VCF file
  const char *vcf_path = argv[optind++];
  htsFile *vcf_file = hts_open(vcf_path, "r");

  bcf_hdr_t *vcfhdr = bcf_hdr_read(vcf_file);
  bcf1_t *v = bcf_init1();
  int s;

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
  FILE *fout = futil_fopen_create(out_path, "w");

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_COVGS);

  // Load reference genome
  ReadBuffer chromsbuf;
  read_buf_alloc(&chromsbuf, 512);
  khash_t(ChromHash) *genome = kh_init(ChromHash);
  seq_reader_load_ref_genome2(ref_buf.b, ref_buf.len, &chromsbuf, genome);

  // Add samples to vcf header
  for(i = 0; i < db_graph.num_of_cols; i++)
    bcf_hdr_add_sample(vcfhdr, db_graph.ginfo[i].sample_name.b);

  // TODO: Load kmers from VCF + ref
  read_t *chrom;
  const char *chrom_name;

  GenoVarBuffer varbuf;
  genovar_buf_alloc(&varbuf, 512);

  while((s = bcf_read(vcf_file, vcfhdr, v)) == 0)
  {
    genovar_buf_reset(&varbuf);
    bcf_unpack(v, BCF_UN_ALL);
    chrom_name = vcfhdr->id[BCF_DT_CTG][v->rid].key;
    chrom = seq_fetch_chrom(genome, chrom_name);
    if(v->n_allele < 2) die("Bad line");
    for(i = 1; i < v->n_allele; i++) {
      genovar_buf_add(&varbuf, (GenoVar){.pos = 0});
      add_new_var(&varbuf.b[varbuf.len-1], v->pos,
                  v->d.allele[0], v->d.allele[i]);
    }

    // fprintf(stderr, "%s %i %u", chrom_name, v->pos+1, v->rlen);
    // fprintf(stderr, " %s", v->d.allele[0]);
    // fprintf(stderr, " %s", v->d.allele[1]);
    // for(i=2; i < v->n_allele; i++)
    //   fprintf(stderr, ",%s", v->d.allele[i]);
    // fprintf(stderr, "\n");
  }

  bcf_destroy(v);
  bcf_hdr_destroy(vcfhdr);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = true,
                              .must_exist_in_edges = NULL,
                              .empty_colours = false};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Seek to the start of VCF file
  // gzseek(vcf_file, SEEK_SET, 0);

  hts_close(vcf_file);
  vcf_file = hts_open(vcf_path, "r");

  // TODO: genotype calls

  status("  saved to: %s\n", out_path);
  
  // gzclose(vcf_file);
  hts_close(vcf_file);
  fclose(fout);

  for(i = 0; i < chromsbuf.len; i++) seq_read_dealloc(&chromsbuf.b[i]);
  read_buf_dealloc(&chromsbuf);
  kh_destroy_ChromHash(genome);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
