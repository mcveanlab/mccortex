#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "gpath_checks.h"
#include "seq_reader.h"

const char uniqkmers_usage[] =
"usage: "CMD" uniqkmers [options] <N>\n"
"\n"
"  Generate <N> random unique kmers not in the graph.\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -o, --out <bub.txt.gz>  Output file [required]\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
//
"  -k, --kmer <K>          Kmer size (required if only giving --seq input)\n"
"  -g, --graph <in.ctx>    Load kmers from the graph file\n"
"  -1, --seq <in.fa>       Load kmers from a sequence file [SAM/BAM/FASTQ etc.]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"paths",        required_argument, NULL, 'p'},
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"kmer",         required_argument, NULL, 'k'},
  {"graph",        required_argument, NULL, 'g'},
  {"seq",          required_argument, NULL, '1'},
  {NULL, 0, NULL, 0}
};

int ctx_uniqkmers(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t i, kmer_size = 0;

  GraphFileReader tmp_gfile;
  GraphFileBuffer gfilebuf;
  gfile_buf_alloc(&gfilebuf, 8);

  seq_file_t *tmp_sfile;
  SeqFilePtrBuffer sfilebuf;
  seq_file_ptr_buf_alloc(&sfilebuf, 16);

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
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'k': cmd_check(!kmer_size,cmd); kmer_size = cmd_kmer_size(cmd, optarg); break;
      case 'g':
        graph_file_reset(&tmp_gfile);
        graph_file_open2(&tmp_gfile, optarg, "r", 0);
        file_filter_flatten(&tmp_gfile.fltr, 0);
        gfile_buf_add(&gfilebuf, tmp_gfile);
        break;
      case '1':
        if((tmp_sfile = seq_open(optarg)) == NULL)
          die("Cannot read --seq file %s", optarg);
        seq_file_ptr_buf_add(&sfilebuf, tmp_sfile);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" uniqkmers -h` for help. Bad option: %s", argv[optind-1]);
      default: die("Bad option: %s", cmd);
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";

  if(optind+1 != argc) cmd_print_usage("Expected <N> for number of kmers only");

  size_t num_uniqkmers = 0;

  if(!parse_entire_size(argv[argc-1], &num_uniqkmers))
    cmd_print_usage("Invalid number of unique kmers: %s", argv[argc-1]);

  if(gfilebuf.len == 0 && !kmer_size)
    die("kmer size not set with -k <K>");
  else if(!kmer_size)
    kmer_size = gfilebuf.data[0].hdr.kmer_size;
  else if(gfilebuf.len > 0 && gfilebuf.data[0].hdr.kmer_size != kmer_size) {
    die("Kmer size mismatches: %zu vs %zu (you do not have to specify -k)",
        (size_t)gfilebuf.data[0].hdr.kmer_size, kmer_size);
  }

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfilebuf.data, gfilebuf.len, NULL, 0, -1);

  bool ctx_uses_stdin = false;
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

  for(i = 0; i < gfilebuf.len; i++) {
    ctx_max_kmers = MAX2(ctx_max_kmers, graph_file_nkmers(&gfilebuf.data[i]));
    ctx_sum_kmers += graph_file_nkmers(&gfilebuf.data[i]);
    ctx_uses_stdin |= file_filter_isstdin(&gfilebuf.data[i].fltr);
  }

  if(ctx_uses_stdin) ctx_sum_kmers = SIZE_MAX;

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  FILE *fout = futil_open_create(out_path, "w");

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, 1, kmers_in_hash, 0);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  for(i = 0; i < gfilebuf.len; i++) {
    graph_load(&gfilebuf.data[i], gprefs, &stats);
    graph_file_close(&gfilebuf.data[i]);
    gprefs.empty_colours = false;
  }

  hash_table_print_stats(&db_graph.ht);

  // Print random kmers
  BinaryKmer bkmer;
  bool found;
  dBNode node;
  char bkmerstr[MAX_KMER_SIZE+1];

  for(i = 0; i < num_uniqkmers; )
  {
    bkmer = binary_kmer_random(kmer_size);
    node = db_graph_find_or_add_node(&db_graph, bkmer, &found);
    if(!found) {
      binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
      fputs(bkmerstr, fout);
      fputc('\n', fout);
      i++;
    }
  }

  char num_kmers_str[100];
  ulong_to_str(num_uniqkmers, num_kmers_str);
  status("  wrote %s kmers to: %s\n", num_kmers_str, futil_outpath_str(out_path));
  fclose(fout);

  gfile_buf_dealloc(&gfilebuf);
  seq_file_ptr_buf_dealloc(&sfilebuf);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
