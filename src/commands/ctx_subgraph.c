#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "seq_reader.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "graph_info.h"
#include "db_node.h"
#include "graphs_load.h"
#include "graph_writer.h"
#include "subgraph.h"
#include "hash_mem.h" // for calculating mem usage

const char subgraph_usage[] =
"usage: "CMD" subgraph [options] <in.ctx>[:cols] [in2.ctx ...]\n"
"\n"
"  Loads graphs (in.ctx) and dumps a graph (out.ctx) that contains all kmers within\n"
"  <dist> edges of kmers in <seeds.fa>.  Maintains number of colours / covgs etc.\n"
"  Loads seed files twice: 1) get seed; 2) extend;  This lowers memory requirement\n"
"  for large (seed) graphs but means seed files cannot be pipes / sockets.\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -f, --force           Overwrite output files\n"
"  -o, --out <out.ctx>   Save output graph file [required]\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
//
"  -N, --ncols <c>       How many colours to load at once [default: 1]\n"
"  -1, --seq <seed.fa>   Read in a seed file [require at least one]\n"
"  -d, --dist <N>        Number of kmers to extend by [default: 0]\n"
// "  -D, --udist <N>       Number of unitigs to extend by [default: 0]\n"
"  -v, --invert          Dump kmers not in subgraph\n"
"  -U, --unitigs         Grab entire runs of kmers that are touched by a read\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"force",        no_argument,       NULL, 'f'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
// command specific
  {"ncols",        required_argument, NULL, 'N'},
  {"seed",         required_argument, NULL, 's'},
  {"seq",          required_argument, NULL, '1'},
  {"dist",         required_argument, NULL, 'd'},
  // {"sdist",        required_argument, NULL, 'D'},
  {"invert",       no_argument,       NULL, 'v'},
  {"unitigs",      no_argument,       NULL, 'U'},
  {NULL, 0, NULL, 0}
};

int ctx_subgraph(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t i, j, use_ncols = 0, dist = 0;
  bool invert = false, grab_unitigs = false;

  seq_file_t *tmp_sfile;
  SeqFilePtrBuffer sfilebuf;
  seq_file_ptr_buf_alloc(&sfilebuf, 16);

  // Arg parsing
  char cmd[100], shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 't': cmd_check(!nthreads,cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'N': cmd_check(!use_ncols,cmd); use_ncols = cmd_uint32_nonzero(cmd, optarg); break;
      case '1':
      case 's':
        if((tmp_sfile = seq_open(optarg)) == NULL)
          die("Cannot read --seq file %s", optarg);
        seq_file_ptr_buf_add(&sfilebuf, tmp_sfile);
        break;
      case 'd': cmd_check(!dist,cmd); dist = cmd_uint32(cmd, optarg); break;
      case 'v': cmd_check(!invert,cmd); invert = true; break;
      case 'S':
      case 'U': cmd_check(!grab_unitigs,cmd); grab_unitigs = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" subgraph -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;
  if(use_ncols == 0) use_ncols = 1;

  if(sfilebuf.len == 0) cmd_print_usage("Require at least one --seq file");
  if(optind >= argc) cmd_print_usage("Require input graph files (.ctx)");

  size_t num_gfiles = argc - optind;
  char **gfile_paths = argv + optind;

  size_t total_cols;

  // Open graph files
  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

  total_cols = graph_files_open(gfile_paths, gfiles, num_gfiles,
                                &ctx_max_kmers, &ctx_sum_kmers);

  if(use_ncols < total_cols && (out_path == NULL || strcmp(out_path,"-")==0))
    cmd_print_usage("Need to use --ncols %zu if output is stdout", total_cols);

  //
  // Decide on memory
  //
  use_ncols = MIN2(use_ncols, total_cols);
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t num_of_fringe_nodes, fringe_mem, total_mem;
  char graph_mem_str[100], fringe_mem_str[100], num_fringe_nodes_str[100];

  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  ((sizeof(Edges) + sizeof(Covg))*use_ncols*8 + 1);

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, &graph_mem);

  graph_mem = hash_table_mem(kmers_in_hash, bits_per_kmer, NULL);
  bytes_to_str(graph_mem, 1, graph_mem_str);

  if(graph_mem >= memargs.mem_to_use)
    die("Not enough memory for graph (requires %s)", graph_mem_str);

  // Fringe nodes
  fringe_mem = memargs.mem_to_use - graph_mem;
  num_of_fringe_nodes = fringe_mem / (sizeof(dBNode) * 2);
  ulong_to_str(num_of_fringe_nodes, num_fringe_nodes_str);
  bytes_to_str(fringe_mem, 1, fringe_mem_str);

  status("[memory] fringe nodes: %s (%s)\n", fringe_mem_str, num_fringe_nodes_str);

  if(dist > 0 && fringe_mem < 1024)
    die("Not enough memory for the graph search (set -m <mem> higher)");

  // Don't need to check, but it prints out memory
  total_mem = graph_mem + fringe_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Open output file
  //

  // Print to stdout unless --out <out> is specified
  if(out_path == NULL) out_path = "-";
  futil_create_output(out_path);

  // Create db_graph
  // multiple colours may be useful later in pulling out multiple colours
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, use_ncols, use_ncols,
                 kmers_in_hash, DBG_ALLOC_EDGES | DBG_ALLOC_COVGS);

  uint8_t *kmer_mask = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

  //
  // Load graphs
  //
  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);

  StrBuf intersect_gname;
  strbuf_alloc(&intersect_gname, 1024);

  if(total_cols > db_graph.num_of_cols) {
    graphs_load_files_flat(gfiles, num_gfiles, gprefs, NULL);
  }
  else {
    for(i = 0; i < num_gfiles; i++)
      graph_load(&gfiles[i], gprefs, NULL);
  }

  // Create header
  for(i = 0; i < num_gfiles; i++) {
    for(j = 0; j < file_filter_num(&gfiles[i].fltr); j++) {
      size_t fromcol = file_filter_fromcol(&gfiles[i].fltr, j);
      graph_info_make_intersect(&gfiles[i].hdr.ginfo[fromcol], &intersect_gname);
    }
  }

  hash_table_print_stats(&db_graph.ht);

  char subgraphstr[] = "subgraph:{";
  strbuf_insert(&intersect_gname, 0, subgraphstr, strlen(subgraphstr));
  strbuf_append_char(&intersect_gname, '}');

  // Load sequence and mark in first pass
  subgraph_from_reads(&db_graph, nthreads, dist,
                      invert, grab_unitigs,
                      fringe_mem, kmer_mask,
                      sfilebuf.b, sfilebuf.len);

  for(i = 0; i < sfilebuf.len; i++) seq_close(sfilebuf.b[i]);
  seq_file_ptr_buf_dealloc(&sfilebuf);

  ctx_free(kmer_mask);
  hash_table_print_stats(&db_graph.ht);

  // Dump nodes that were flagged
  Edges *intersect_edges = NULL;
  bool kmers_loaded = true;
  bool colours_loaded = (total_cols <= db_graph.num_of_cols);

  if(!colours_loaded)
  {
    // Need to reload graph colours - therefore construct edge intersection set
    intersect_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));
    for(i = 0; i < db_graph.ht.capacity; i++)
      intersect_edges[i] = db_node_get_edges_union(&db_graph, i);
  }

  graph_writer_merge_mkhdr(out_path, gfiles, num_gfiles,
                          kmers_loaded, colours_loaded,
                          intersect_edges, intersect_gname.b,
                          &db_graph);

  ctx_free(intersect_edges);
  strbuf_dealloc(&intersect_gname);
  for(i = 0; i < num_gfiles; i++) graph_file_close(&gfiles[i]);
  ctx_free(gfiles);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
