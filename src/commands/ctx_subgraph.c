#include "global.h"

#include "seq_file.h"

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

const char subgraph_usage[] =
"usage: "CMD" subgraph [options] <out.ctx> <dist> <in.ctx>[:cols] [in2.ctx ...]\n"
"\n"
"  Loads graphs (in.ctx) and dumps a graph (out.ctx) that contains all kmers within\n"
"  <dist> edges of kmers in <seeds.fa>.  Maintains number of colours / covgs etc.\n"
"  Loads seed files twice: 1) get seed; 2) extend;  This lowers memory requirement\n"
"  for large (seed) graphs but means seed files cannot be pipes / sockets.\n"
"\n"
"  Options:\n"
"    -m <mem>          Memory to use  <required>\n"
"    -n <kmers>        Hash size\n"
"    --seq <seed.fa>   Read in a seed file\n"
"    --invert          Dump kmers not in subgraph\n"
"    --ncols <n>       Number of samples in memory at once (speedup)\n";

int ctx_subgraph(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that we have at least 4 args

  seq_file_t *seed_files[argc/2];
  size_t num_seed_files = 0;
  bool invert = false;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(!strcasecmp(argv[argi], "--seq") | !strcasecmp(argv[argi], "--seed"))
    {
      if(argi+1 == argc)
        cmd_print_usage("%s <seed.fa> requires an argument", argv[argi]);
      seed_files[num_seed_files] = seq_open(argv[argi+1]);
      if(seed_files[num_seed_files] == NULL)
        die("Cannot read %s file: %s", argv[argi], argv[argi+1]);
      argi++; num_seed_files++;
    }
    else if(strcasecmp(argv[argi], "--invert") == 0) invert = true;
    else cmd_print_usage("Unknown option: %s", argv[argi]);
  }

  const char *out_path = argv[argi], *diststr = argv[argi+1];
  uint32_t dist;

  if(futil_file_exists(out_path))
    die("Output file already exists: %s", out_path);

  if(!parse_entire_uint(diststr, &dist))
    cmd_print_usage("Invalid <dist> value, must be int >= 0: %s", diststr);

  int num_gfiles_int = argc - 2*(int)num_seed_files - 2;
  if(num_gfiles_int <= 0)
    cmd_print_usage("Please specify input graph files (.ctx)");

  size_t i, j, col, num_gfiles = (size_t)num_gfiles_int, total_cols = 0;
  char **paths = argv + 2*num_seed_files + 2;

  // Open graph files
  uint64_t max_num_kmers = 0;
  GraphFileReader gfiles[num_gfiles];

  for(i = 0; i < num_gfiles; i++)
  {
    gfiles[i] = INIT_GRAPH_READER;
    graph_file_open(&gfiles[i], paths[i], true);

    if(gfiles[0].hdr.kmer_size != gfiles[i].hdr.kmer_size) {
      die("Graph kmer-sizes do not match [%u vs %u; %s; %s]\n",
          gfiles[0].hdr.kmer_size, gfiles[i].hdr.kmer_size,
          gfiles[0].fltr.file_path.buff, gfiles[i].fltr.file_path.buff);
    }

    size_t offset = total_cols;
    total_cols += graph_file_usedcols(&gfiles[i]);
    file_filter_update_intocol(&gfiles[i].fltr, gfiles[i].fltr.intocol + offset);

    max_num_kmers = MAX2(gfiles[i].num_of_kmers, max_num_kmers);
  }

  //
  // Decide on memory
  //
  const size_t use_ncols = MIN2(args->use_ncols, total_cols);
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t num_of_fringe_nodes, fringe_mem, total_mem;
  char graph_mem_str[100], fringe_mem_str[100], num_fringe_nodes_str[100];

  bits_per_kmer = ((sizeof(Edges) + sizeof(Covg))*use_ncols*8 + 1);
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, max_num_kmers,
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
  cmd_check_mem_limit(args, total_mem);

  // Check output directory
  if(!futil_is_file_writable(out_path))
    die("Cannot write to output file: %s", out_path);

  // Create db_graph
  // multiple colours may be useful later in pulling out multiple colours
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, use_ncols, use_ncols, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity*use_ncols, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity*use_ncols, sizeof(Covg));

  size_t num_words64 = roundup_bits2words64(db_graph.ht.capacity);
  uint64_t *kmer_mask = calloc2(num_words64, sizeof(uint64_t));

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
  subgraph_from_reads(&db_graph, dist, invert, fringe_mem, kmer_mask,
                      seed_files, num_seed_files);

  for(i = 0; i < num_seed_files; i++) seq_close(seed_files[i]);

  free(kmer_mask);
  hash_table_print_stats(&db_graph.ht);

  // Dump nodes that were flagged
  Edges *intersect_edges = NULL;
  bool kmers_loaded = true;
  bool colours_loaded = (total_cols <= db_graph.num_of_cols);

  if(!colours_loaded)
  {
    // Need to reload graph colours - therefore construct edge intersection set
    intersect_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
    for(i = 0; i < db_graph.ht.capacity; i++)
      intersect_edges[i] = db_node_get_edges_union(&db_graph, i);
  }

  graph_files_merge_mkhdr(out_path, gfiles, num_gfiles,
                          kmers_loaded, colours_loaded,
                          intersect_edges, intersect_gname.buff,
                          &db_graph);

  if(intersect_edges != NULL) free(intersect_edges);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);

  for(i = 0; i < num_gfiles; i++) graph_file_dealloc(&gfiles[i]);

  strbuf_dealloc(&intersect_gname);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
