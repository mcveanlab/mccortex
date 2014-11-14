#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "assemble_contigs.h"
#include "seq_reader.h"
#include "graph_format.h"
#include "gpath_reader.h"
#include "gpath_checks.h"

const char contigs_usage[] =
"usage: "CMD" contigs [options] <input.ctx>\n"
"\n"
"  Assemble contigs from the graph, print statistics\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -f, --force           Overwrite output files\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <N>      Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -o, --out <out.fa>    Print contigs in FASTA [default: don't print]\n"
"  -c, --colour <c>      Pull out contigs from the given colour [default: 0]\n"
"  -p, --paths <in.ctp>  Load path file (can specify multiple times)\n"
"  -N, --ncontigs <N>    Pull out <N> contigs from random kmers [default: 0, no limit]\n"
"  -s, --seed <in.fa>    Use seed kmers from a file. Reads must be of kmer length\n"
"  -r, --reseed          Sample seed kmers with replacement\n"
"  -R, --no-reseed       Do not use a seed kmer if it is used in a contig [default]\n"
"  -P, --use-seed-paths  Use unused paths to seed contigs [default: off]\n"
"  -G, --genome <G>      Genome size in bases\n"
"  -S, --confid-csv <save.csv> Save confidence table to <save.csv>\n"
"  -M, --no-missing-check      Do not use the missing information check\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"paths",        required_argument, NULL, 'p'},
// command specific
  {"seed",         required_argument, NULL, 's'},
  {"seq",          required_argument, NULL, '1'},
  {"reseed",       no_argument,       NULL, 'r'},
  {"no-reseed",    no_argument,       NULL, 'R'},
  {"use-seed-paths",no_argument,      NULL, 'P'},
  {"ncontigs",     required_argument, NULL, 'N'},
  {"colour",       required_argument, NULL, 'c'},
  {"color",        required_argument, NULL, 'c'},
  {"genome",       required_argument, NULL, 'G'},
  {"confid-csv",   required_argument, NULL, 'S'},
  {"no-missing-check", no_argument,   NULL, 'M'},
  {NULL, 0, NULL, 0}
};

int ctx_contigs(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t i, contig_limit = 0, colour = 0;
  bool cmd_reseed = false, cmd_no_reseed = false; // -r, -R
  const char *conf_table_path = NULL; // save confidence table to here
  bool use_missing_info_check = true, seed_with_unused_paths = false;

  // Read length and expected depth for calculating confidences
  size_t genome_size = 0;

  seq_file_t *tmp_seed_file = NULL;
  SeqFilePtrBuffer seed_buf;
  seq_file_ptr_buf_alloc(&seed_buf, 16);

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

  // Arg parsing
  char cmd[100], shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!out_path,cmd); out_path = optarg; break;
      case 't': cmd_check(!nthreads,cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_add(&gpfiles, tmp_gpfile);
        break;
      case '1':
      case 's': // --seed <in.fa>
        if((tmp_seed_file = seq_open(optarg)) == NULL)
          die("Cannot read --seed file: %s", optarg);
        seq_file_ptr_buf_add(&seed_buf, tmp_seed_file);
        break;
      case 'r': cmd_check(!cmd_reseed,cmd); cmd_reseed = true; break;
      case 'R': cmd_check(!cmd_no_reseed,cmd); cmd_no_reseed = true; break;
      case 'N':
        cmd_check(!contig_limit,cmd);
        contig_limit = cmd_uint32_nonzero(cmd, optarg);
        break;
      case 'c': cmd_check(!colour,cmd); colour = cmd_uint32(cmd, optarg); break;
      case 'G': cmd_check(!genome_size,cmd); genome_size = cmd_bases(cmd, optarg); break;
      case 'S': cmd_check(!conf_table_path,cmd); conf_table_path = optarg; break;
      case 'M': cmd_check(use_missing_info_check,cmd); use_missing_info_check = false; break;
      case 'P': cmd_check(!seed_with_unused_paths,cmd); seed_with_unused_paths = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" contigs -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(cmd_no_reseed && cmd_reseed)
    cmd_print_usage("Cannot specify both -r and -R");

  bool sample_with_replacement = cmd_reseed;

  // Defaults
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;

  if(!seed_buf.len && !contig_limit && sample_with_replacement) {
    cmd_print_usage("Please specify one or more of: "
                    "--no-reseed | --ncontigs | --seed <in.fa>");
  }

  if(optind+1 != argc)
    cmd_print_usage("Too %s arguments", optind == argc ? "few" : "many");

  char *ctx_path = argv[optind];

  //
  // Open Graph file
  //
  GraphFileReader gfile;
  memset(&gfile, 0, sizeof(GraphFileReader));
  graph_file_open(&gfile, ctx_path);

  // Update colours in graph file - sample in 0, all others in 1
  // never need more than two colours
  size_t ncols = gpath_load_sample_pop(&gfile, gpfiles.data, gpfiles.len, colour);

  // Check for compatibility between graph files and path files
  // pop_colour is colour 1
  graphs_gpaths_compatible(&gfile, 1, gpfiles.data, gpfiles.len, 1);

  if(!genome_size)
  {
    char nk_str[50];
    if(gfile.num_of_kmers < 0) die("Please pass --genome <G> if streaming");
    genome_size = gfile.num_of_kmers;
    ulong_to_str(genome_size, nk_str);
    status("Taking number of kmers as genome size: %s", nk_str);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;

  // 1 bit needed per kmer if we need to keep track of kmer usage
  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 + sizeof(GPath*)*8 +
                  ncols + !sample_with_replacement;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        gfile.num_of_kmers, gfile.num_of_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.data, gpfiles.len, ncols, rem_mem, false);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  // Total memory
  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  // Load contig hist distribution from ctp files
  ZeroSizeBuffer contig_hist;
  memset(&contig_hist, 0, sizeof(contig_hist));

  for(i = 0; i < gpfiles.len; i++) {
    gpath_reader_load_contig_hist(gpfiles.data[i].json,
                                  gpfiles.data[i].fltr.path.b,
                                  file_filter_fromcol(&gpfiles.data[i].fltr, 0),
                                  &contig_hist);
  }

  // Calculate confidences, only for one colour
  ContigConfidenceTable conf_table;
  conf_table_alloc(&conf_table, 1);
  conf_table_update_hist(&conf_table, 0, genome_size,
                         contig_hist.data, contig_hist.len);

  if(conf_table_path != NULL) {
    conf_table_save(&conf_table, conf_table_path);
  }

  zsize_buf_dealloc(&contig_hist);

  //
  // Output file if printing
  //
  FILE *fout = out_path ? futil_open_create(out_path, "w") : NULL;

  // Allocate
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfile.hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.data, gpfiles.len, path_mem,
                             false, &db_graph);

  uint8_t *visited = NULL;

  if(!sample_with_replacement)
    visited = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

  // Load graph
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = true};

  graph_load(&gfile, gprefs, &stats);
  graph_file_close(&gfile);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  for(i = 0; i < gpfiles.len; i++) {
    gpath_reader_load(&gpfiles.data[i], GPATH_DIE_MISSING_KMERS, &db_graph);
    gpath_reader_close(&gpfiles.data[i]);
  }
  gpfile_buf_dealloc(&gpfiles);

  AssembleContigStats assem_stats;
  assemble_contigs_stats_init(&assem_stats);

  assemble_contigs(nthreads, seed_buf.data, seed_buf.len,
                   contig_limit, visited,
                   use_missing_info_check, seed_with_unused_paths,
                   fout, out_path, &assem_stats, &conf_table,
                   &db_graph, 0); // Sample always loaded into colour zero

  if(fout && fout != stdout) fclose(fout);

  assemble_contigs_stats_print(&assem_stats);
  assemble_contigs_stats_destroy(&assem_stats);

  conf_table_dealloc(&conf_table);

  for(i = 0; i < seed_buf.len; i++)
    seq_close(seed_buf.data[i]);

  seq_file_ptr_buf_dealloc(&seed_buf);

  ctx_free(visited);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
