#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "graph_info.h"
#include "db_node.h"
#include "graph_format.h"
#include "subgraph.h"

#include "seq_file.h"

const char subgraph_usage[] =
"usage: "CMD" subgraph [options] <in.ctx>[:cols] [in2.ctx ...]\n"
"\n"
"  Loads graphs (in.ctx) and dumps a graph (out.ctx) that contains all kmers within\n"
"  <dist> edges of kmers in <seeds.fa>.  Maintains number of colours / covgs etc.\n"
"  Loads seed files twice: 1) get seed; 2) extend;  This lowers memory requirement\n"
"  for large (seed) graphs but means seed files cannot be pipes / sockets.\n"
"\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -s, --ncols <c>       How many colours to load at once [default: 1]\n"
"  -o, --out <out.ctx>   Save output graph file [required]\n"
"  -1, --seq <seed.fa>   Read in a seed file [require at least one]\n"
"  -d, --dist <N>        Number of kmers to extend by [default: 0]\n"
// "  -D, --sdist <N>       Number of supernodes to extend by [default: 0]\n"
"  -v, --invert          Dump kmers not in subgraph\n"
"  -u, --supernodes      Grab entire runs of kmers that are touched by a read\n"
"\n";

int ctx_subgraph(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that we have at least 4 args

  seq_file_t *seed_files[argc/2];
  size_t num_seed_files = 0, dist = 0;
  bool invert = false, grab_supernodes = false;

  size_t num_threads = args->max_work_threads;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++)
  {
    if(!strcmp(argv[argi], "--seq") | !strcmp(argv[argi], "--seed"))
    {
      if(argi+1 == argc)
        cmd_print_usage("%s <seed.fa> requires an argument", argv[argi]);
      seed_files[num_seed_files] = seq_open(argv[argi+1]);
      if(seed_files[num_seed_files] == NULL)
        die("Cannot read %s file: %s", argv[argi], argv[argi+1]);
      argi++; num_seed_files++;
    }
    else if(!strcmp(argv[argi], "--dist") || !strcmp(argv[argi], "-d")) {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &dist)) {
        cmd_print_usage("%s <N> requires an integer argument >= 0", argv[argi]);
      }
      argi++; // we took an argument
    }
    else if(!strcmp(argv[argi], "--invert") || !strcmp(argv[argi], "-v"))
      invert = true;
    else if(!strcmp(argv[argi], "--supernodes") || !strcmp(argv[argi], "-u"))
      grab_supernodes = true;
    else cmd_print_usage("Unknown option: %s", argv[argi]);
  }

  if(argi >= argc)
    cmd_print_usage("Please specify at least one input graph file (.ctx)");

  size_t num_gfiles = argc - argi;
  char **gfile_paths = argv + argi;

  size_t i, j, col, total_cols;

  // Open graph files
  GraphFileReader gfiles[num_gfiles];
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

  total_cols = graph_files_open(gfile_paths, gfiles, num_gfiles,
                                &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Decide on memory
  //
  const size_t use_ncols = MIN2(args->use_ncols, total_cols);
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t num_of_fringe_nodes, fringe_mem, total_mem;
  char graph_mem_str[100], fringe_mem_str[100], num_fringe_nodes_str[100];

  bits_per_kmer = ((sizeof(Edges) + sizeof(Covg))*use_ncols*8 + 1);
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, NULL);

  graph_mem = hash_table_mem(kmers_in_hash,bits_per_kmer,NULL);
  bytes_to_str(graph_mem, 1, graph_mem_str);

  if(graph_mem >= args->mem_to_use)
    die("Not enough memory for graph (requires %s)", graph_mem_str);

  // Fringe nodes
  fringe_mem = args->mem_to_use - graph_mem;
  num_of_fringe_nodes = fringe_mem / (sizeof(dBNode) * 2);
  ulong_to_str(num_of_fringe_nodes, num_fringe_nodes_str);
  bytes_to_str(fringe_mem, 1, fringe_mem_str);

  status("[memory] fringe nodes: %s (%s)\n", fringe_mem_str, num_fringe_nodes_str);

  if(dist > 0 && fringe_mem < 1024)
    die("Not enough memory for the graph search (set -m <mem> higher)");

  // Don't need to check, but it prints out memory
  total_mem = graph_mem + fringe_mem;
  cmd_check_mem_limit(args->mem_to_use, total_mem);

  //
  // Open output file
  //

  // Print to stdout unless --out <out> is specified
  const char *out_path = args->output_file_set ? args->output_file : "-";

  if(args->output_file_set) {
    if(futil_file_exists(out_path)) die("File already exists: %s", out_path);
    else if(!futil_is_file_writable(out_path)) die("Cannot write: %s", out_path);
  }

  // Create db_graph
  // multiple colours may be useful later in pulling out multiple colours
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, use_ncols, use_ncols, kmers_in_hash);
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity*use_ncols, sizeof(Edges));
  db_graph.col_covgs = ctx_calloc(db_graph.ht.capacity*use_ncols, sizeof(Covg));

  uint8_t *kmer_mask = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  //
  // Load graphs
  //
  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  StrBuf intersect_gname;
  strbuf_alloc(&intersect_gname, 1024);

  size_t tmpinto; bool tmpflatten;
  for(i = 0; i < num_gfiles; i++) {
    tmpinto = gfiles[i].fltr.intocol;
    tmpflatten = gfiles[i].fltr.flatten;

    if(total_cols > db_graph.num_of_cols) {
      file_filter_update_intocol(&gfiles[i].fltr, 0);
      gfiles[i].fltr.flatten = true;
    }

    graph_load(&gfiles[i], gprefs, &stats);
    file_filter_update_intocol(&gfiles[i].fltr, tmpinto);
    gfiles[i].fltr.flatten = tmpflatten;

    for(j = 0; j < gfiles[i].fltr.ncols; j++) {
      col = gfiles[i].fltr.cols[j];
      graph_info_make_intersect(&gfiles[i].hdr.ginfo[col], &intersect_gname);
    }
  }

  hash_table_print_stats(&db_graph.ht);

  char subgraphstr[] = "subgraph:{";
  strbuf_insert(&intersect_gname, 0, subgraphstr, strlen(subgraphstr));
  strbuf_append_char(&intersect_gname, '}');

  // Load sequence and mark in first pass
  subgraph_from_reads(&db_graph, num_threads, dist,
                      invert, grab_supernodes,
                      fringe_mem, kmer_mask,
                      seed_files, num_seed_files);

  for(i = 0; i < num_seed_files; i++) seq_close(seed_files[i]);

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

  graph_files_merge_mkhdr(out_path, gfiles, num_gfiles,
                          kmers_loaded, colours_loaded,
                          intersect_edges, intersect_gname.buff,
                          &db_graph);

  ctx_free(intersect_edges);
  strbuf_dealloc(&intersect_gname);
  for(i = 0; i < num_gfiles; i++) graph_file_close(&gfiles[i]);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
