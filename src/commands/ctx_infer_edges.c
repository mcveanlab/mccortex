#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "graph_file_reader.h"
#include "loading_stats.h"
#include "seq_reader.h"
#include "infer_edges.h"

const char inferedges_usage[] =
"usage: "CMD" inferedges [options] <pop.ctx>\n"
"\n"
"  Infer edges in a population graph.  By default adds all missing edges (--all).\n"
"  It is important that you run this step before doing read threading.\n"
"\n"
"  -h, --help            This help message\n"
"  -o, --out <out.ctp>   Save output file\n"
"  -m, --memory <mem>    Memory to use (e.g. 1M, 20GB)\n"
"  -n, --nkmers <N>      Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -P, --pop             Add edges that are in the union only\n"
"  -A, --all             Add all edges [default]\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
// command specific
  {"pop",          no_argument,       NULL, 'P'},
  {"all",          no_argument,       NULL, 'A'},
  {NULL, 0, NULL, 0}
};

// If two kmers are in a sample and the population has an edges between them,
// Add edge to sample

// Using file so can call fseek and don't need to load whole graph
static size_t inferedges_on_file(const dBGraph *db_graph, bool add_all_edges,
                                 GraphFileReader *file, FILE *fout)
{
  ctx_assert(file->hdr.num_of_cols == db_graph->num_of_cols);
  ctx_assert2(file->fltr.fh != stdin, "Use inferedges_on_stream() instead");
  const size_t ncols = db_graph->num_of_cols;

  status("[inferedges] Processing file: %s", file->fltr.file_path.buff);

  if(fout != NULL) {
    // Print header
    graph_write_header(fout, &file->hdr);
  }

  // Read the input file again
  if(fseek(file->fltr.fh, file->hdr_size, SEEK_SET) != 0)
    die("fseek failed: %s", strerror(errno));

  BinaryKmer bkmer;
  Edges edges[ncols];
  Covg covgs[ncols];

  size_t num_nodes_modified = 0;
  bool updated;
  long edges_len = sizeof(Edges) * ncols;

  while(graph_file_read(file, &bkmer, covgs, edges))
  {
    updated = (add_all_edges ? infer_all_edges(bkmer, edges, covgs, db_graph)
                             : infer_pop_edges(bkmer, edges, covgs, db_graph));

    if(fout != NULL) {
      graph_write_kmer(fout, file->hdr.num_of_bitfields, file->hdr.num_of_cols,
                       bkmer, covgs, edges);
    }
    else if(updated) {
      if(fseek(file->fltr.fh, -edges_len, SEEK_CUR) != 0 ||
         fwrite(edges, sizeof(Edges), ncols, file->fltr.fh) != ncols)
      {
        die("fseek/fwrite error: %s", file->fltr.file_path.buff);
      }
    }

    num_nodes_modified += updated;
  }

  return num_nodes_modified;
}


int ctx_infer_edges(int argc, char **argv)
{
  size_t num_of_threads = DEFAULT_NTHREADS;
  struct MemArgs memargs = MEM_ARGS_INIT;
  char *out_ctx_path = NULL;
  bool add_pop_edges = false, add_all_edges = false;

  // Arg parsing
  char cmd[100];
  char shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': if(out_ctx_path){cmd_print_usage(NULL);} out_ctx_path = optarg; break;
      case 't': num_of_threads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'A': add_all_edges = true; break;
      case 'P': add_pop_edges = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" inferedges -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Default to adding all edges
  if(!add_pop_edges && !add_all_edges) add_all_edges = true;

  // Can only specify one of --pop --all
  if(add_pop_edges && add_all_edges)
    cmd_print_usage("Please specify only one of --all --pop");

  // Check that optind+1 == argc
  if(optind+1 > argc)
    cmd_print_usage("Expected exactly one graph file");
  else if(optind+1 < argc)
    cmd_print_usage("Expected only one graph file. What is this: '%s'", argv[optind]);

  //
  // Open graph file
  //
  char *graph_path = argv[optind];
  status("Reading graph: %s", graph_path);

  GraphFileReader file = INIT_GRAPH_READER;
  // Mode r+ means open (not create) for update (read & write)
  graph_file_open2(&file, graph_path, true, "r+");
  bool reading_stream = (file.fltr.fh == stdin);

  FILE *fout = NULL;

  if(out_ctx_path || reading_stream)
    fout = futil_open_output(out_ctx_path ? out_ctx_path : "-");

  if(!file.fltr.nofilter)
    cmd_print_usage("Inferedges with filter not implemented - sorry");

  if(!futil_is_file_writable(graph_path))
    cmd_print_usage("Cannot write to file: %s", graph_path);

  ctx_assert(file.hdr.num_of_cols == file.fltr.ncols);

  // Print output status
  if(fout == stdout) status("Writing to STDOUT");
  else if(fout != NULL) status("Writing to: %s", out_ctx_path);
  else status("Editing file in place: %s", graph_path);

  status("Inferring all missing %sedges", add_pop_edges ? "population " : "");

  //
  // Decide on memory
  //
  size_t kmers_in_hash, graph_mem, bits_per_kmer;

  // reading file: one bit per kmer per colour: for 'in colour'
  // reading stream: 9 bits per kmer per colour: Edges + one bit for 'in colour'.
  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  file.hdr.num_of_cols * (sizeof(Edges)*8*reading_stream + 1);

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        file.num_of_kmers, file.num_of_kmers,
                                        memargs.mem_to_use_set, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, file.hdr.kmer_size,
                 file.hdr.num_of_cols, file.hdr.num_of_cols*reading_stream,
                 kmers_in_hash);

  // In colour
  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);
  db_graph.node_in_cols = ctx_calloc(bytes_per_col*file.hdr.num_of_cols, 1);

  if(reading_stream)
    db_graph.col_edges = ctx_calloc(file.hdr.num_of_cols*db_graph.ht.capacity, 1);

  LoadingStats stats = LOAD_STATS_INIT_MACRO;
  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = false};

  // We need to load the graph for both --pop and --all since we need to check
  // if the next kmer is in each of the colours
  graph_load(&file, gprefs, &stats);

  if(add_pop_edges) status("Inferring edges from population...\n");
  else status("Inferring all missing edges...\n");

  size_t num_nodes_modified;

  if(reading_stream)
  {
    num_nodes_modified = infer_edges(num_of_threads, add_all_edges, &db_graph);
    graph_write_header(fout, &file.hdr);
    graph_write_all_kmers(fout, &db_graph);
  }
  else {
    num_nodes_modified = inferedges_on_file(&db_graph, add_all_edges,
                                            &file, fout);
  }

  if(fout != NULL && fout != stdout) fclose(fout);

  char modified_str[100], kmers_str[100];
  ulong_to_str(num_nodes_modified, modified_str);
  ulong_to_str(db_graph.ht.num_kmers, kmers_str);

  double modified_rate = 0;
  if(db_graph.ht.num_kmers)
    modified_rate = (100.0 * num_nodes_modified) / db_graph.ht.num_kmers;

  status("%s of %s (%.2f%%) nodes modified\n",
         modified_str, kmers_str, modified_rate);

  graph_file_close(&file);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
