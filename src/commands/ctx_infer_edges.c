#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "loading_stats.h"
#include "seq_reader.h"
#include "infer_edges.h"

const char inferedges_usage[] =
"usage: "CMD" inferedges [options] <pop.ctx>\n"
"\n"
"  Infer edges in a population graph.  One of -P, -A required.\n"
"\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use\n"
"  -o, --out <out.ctx>   Output file [default is to edit input file]\n"
"  -P, --pop             Add edges that are in the union only\n"
"  -A, --all             Add all edges [default]\n"
"\n";

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


int ctx_infer_edges(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that: 2<= argc <=3

  bool add_pop_edges = false, add_all_edges = false;
  while(argc > 0 && argv[0][0] == '-' && argv[0][1]) {
    if(!strcmp(argv[0],"--all") || !strcmp(argv[0],"-A")) {
      argc--; argv++; add_all_edges = true;
    }
    else if(!strcmp(argv[0],"--pop") || !strcmp(argv[0],"-P")) {
      argc--; argv++; add_pop_edges = true;
    }
  }

  // Default to adding all edges
  if(!add_pop_edges && !add_all_edges) add_all_edges = true;

  // Can only specify one of --pop --all
  if(add_pop_edges && add_all_edges)
    cmd_print_usage("Please specify only one of --all --pop");

  if(argc != 1) cmd_print_usage(NULL);

  char *path = argv[0];
  dBGraph db_graph;
  GraphFileReader file = INIT_GRAPH_READER;
  // Mode r+ means open (not create) for update (read & write)
  graph_file_open2(&file, path, true, "r+");
  bool reading_stream = (file.fltr.fh == stdin);
  size_t num_threads = args->max_work_threads;

  FILE *fout = NULL;
  if(args->output_file_set) {
    if(strcmp(args->output_file,"-") == 0) fout = stdout;
    else if(futil_file_exists(args->output_file))
      die("Output file already exists: %s", args->output_file);
    else if((fout = fopen(args->output_file,"w")) == NULL)
      die("Cannot open output file: %s", args->output_file);
    setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);
  }
  else if(reading_stream)
    fout = stdout;
  else
    status("Editing file in place: %s", path);

  if(!file.fltr.nofilter)
    cmd_print_usage("Inferedges with filter not implemented - sorry");

  if(!futil_is_file_writable(path))
    cmd_print_usage("Cannot write to file: %s", path);

  ctx_assert(file.hdr.num_of_cols == file.fltr.ncols);

  // Print output status
  if(fout == stdout) status("Writing to STDOUT");
  else if(fout != NULL) status("Writing to: %s", args->output_file);
  else status("Editing file in place: %s", path);

  //
  // Decide on memory
  //
  size_t kmers_in_hash, graph_mem, extra_bits_per_kmer;

  // reading file: one bit per kmer per colour: for 'in colour'
  // reading stream: 9 bits per kmer per colour: Edges + one bit for 'in colour'.
  extra_bits_per_kmer = file.hdr.num_of_cols * (sizeof(Edges)*8*reading_stream + 1);

  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        file.num_of_kmers, file.num_of_kmers,
                                        args->mem_to_use_set, &graph_mem);

  cmd_check_mem_limit(args->mem_to_use, graph_mem);

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
    num_nodes_modified = infer_edges(num_threads, add_all_edges, &db_graph);
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

  if(reading_stream) ctx_free(db_graph.col_edges);
  ctx_free(db_graph.node_in_cols);
  db_graph_dealloc(&db_graph);

  graph_file_close(&file);

  return EXIT_SUCCESS;
}
