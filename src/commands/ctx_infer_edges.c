#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graphs_load.h"
#include "graph_writer.h"
#include "seq_reader.h"
#include "infer_edges.h"

// Memory mapped files used in inferedges_on_mmap()
#include <sys/mman.h>

const char inferedges_usage[] =
"usage: "CMD" inferedges [options] <pop.ctx>\n"
"\n"
"  Infer edges adds edges between all kmers that share k-1 bases.\n"
"  By default adds all missing edges (--all). To add only edges that exist in\n"
"  at least one other sample in the population, use --pop.\n"
"  It is important that you run this step before doing read threading.\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -f, --force           Overwrite output files\n"
"  -o, --out <out.ctx>   Save output graph file\n"
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
  {"force",        no_argument,       NULL, 'f'},
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
static size_t inferedges_on_mmap(const dBGraph *db_graph, bool add_all_edges,
                                 GraphFileReader *file)
{
  // ctx_assert2(file->strm.b == NULL, "File should not be buffered");
  ctx_assert(db_graph->num_of_cols == file->hdr.num_of_cols);
  ctx_assert(file_filter_from_direct(&file->fltr));
  ctx_assert2(!isatty(fileno(file->fh)), "Use inferedges_on_stream() instead");
  ctx_assert(file->num_of_kmers >= 0);
  ctx_assert(file->file_size >= 0);

  status("[inferedges] Processing mmap file: %s [hdr: %zu bytes file: %zu bytes]",
         file_filter_path(&file->fltr),
         (size_t)file->hdr_size, (size_t)file->file_size);

  if(graph_file_fseek(file, 0, SEEK_SET) != 0)
    die("fseek failed: %s", strerror(errno));

  // Open memory mapped file
  void *mmap_ptr = mmap(NULL, file->file_size, PROT_WRITE, MAP_SHARED,
                        fileno(file->fh), 0);

  if(mmap_ptr == MAP_FAILED)
    die("Cannot memory map file: %s [%s]", file->fltr.path.b, strerror(errno));

  const size_t ncols = file->hdr.num_of_cols;
  BinaryKmer bkmer;
  Edges edges[ncols];
  Covg covgs[ncols];

  bool updated;
  size_t i, num_kmers = file->num_of_kmers, num_kmers_edited = 0;
  size_t filekmersize = sizeof(BinaryKmer) + (sizeof(Edges)+sizeof(Covg)) * ncols;

  char *ptr = (char*)mmap_ptr + file->hdr_size;

  for(i = 0; i < num_kmers; i++, ptr += filekmersize)
  {
    char *fh_covgs = ptr      + sizeof(BinaryKmer);
    char *fh_edges = fh_covgs + sizeof(Covg)*ncols;

    memcpy(bkmer.b, ptr,      sizeof(BinaryKmer));
    memcpy(covgs,   fh_covgs, ncols * sizeof(Covg));
    memcpy(edges,   fh_edges, ncols * sizeof(Edges));

    updated = infer_kmer_edges(bkmer, !add_all_edges, edges, covgs, db_graph);

    if(updated) {
      memcpy(fh_covgs, covgs, ncols * sizeof(Covg));
      memcpy(fh_edges, edges, ncols * sizeof(Edges));
      num_kmers_edited++;
    }
  }

  if(munmap(mmap_ptr, file->file_size) == -1)
    die("Cannot release mmap file: %s [%s]", file->fltr.path.b, strerror(errno));

  return num_kmers_edited;
}


// Using file so can call fseek and don't need to load whole graph
static size_t inferedges_on_file(const dBGraph *db_graph, bool add_all_edges,
                                 GraphFileReader *file, FILE *fout)
{
  // ctx_assert(db_graph->num_of_cols == file->hdr.num_of_cols);
  ctx_assert(file_filter_from_direct(&file->fltr));
  ctx_assert2(!isatty(fileno(file->fh)), "Use inferedges_on_stream() instead");
  ctx_assert(fout != NULL);
  ctx_assert(fileno(file->fh) != fileno(fout));

  status("[inferedges] Processing file: %s", file_filter_path(&file->fltr));

  // Print header
  graph_write_header(fout, &file->hdr);

  // Read the input file again
  if(graph_file_fseek(file, file->hdr_size, SEEK_SET) != 0)
    die("graph_file_fseek failed: %s", strerror(errno));

  const size_t file_ncols = file->hdr.num_of_cols;
  BinaryKmer bkmer;
  Edges edges[file_ncols];
  Covg covgs[file_ncols];

  size_t num_kmers_edited = 0;
  bool updated;

  while(graph_file_read_reset(file, &bkmer, covgs, edges))
  {
    updated = infer_kmer_edges(bkmer, !add_all_edges, edges, covgs, db_graph);
    graph_write_kmer(fout, file_ncols, bkmer, covgs, edges);

    num_kmers_edited += updated;
  }

  return num_kmers_edited;
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
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!out_ctx_path,cmd); out_ctx_path = optarg; break;
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
  status("Reading graph: %s", (!strcmp(graph_path,"-") ? "STDIN" : graph_path));

  if(strchr(graph_path,':') != NULL)
    cmd_print_usage("Cannot use ':' in input graph for `"CMD" inferedges`");

  GraphFileReader file;
  memset(&file, 0, sizeof(file));

  file_filter_open(&file.fltr, graph_path);

  // Use stat to detect if we are reading from a stream
  struct stat st;
  bool reading_stream = (stat(file.fltr.path.b, &st) != 0);
  bool editing_file = !(out_ctx_path || reading_stream);

  // Mode r+ means open (not create) for update (read & write)
  bool use_buf = true;
  graph_file_open2(&file, graph_path, reading_stream ? "r" : "r+", use_buf, 0);

  if(!file_filter_from_direct(&file.fltr))
    cmd_print_usage("Inferedges with filter not implemented - sorry");

  FILE *fout = NULL;

  // Editing input file or writing a new file
  if(!editing_file)
    fout = futil_fopen_create(out_ctx_path ? out_ctx_path : "-", "w");

  // Print output status
  if(fout == stdout) status("Writing to STDOUT");
  else if(fout != NULL) status("Writing to: %s", out_ctx_path);
  else status("Editing file in place: %s", graph_path);

  status("Inferring all missing %sedges", add_pop_edges ? "population " : "");

  //
  // Decide on memory
  //
  const size_t ncols = file.hdr.num_of_cols;
  size_t kmers_in_hash, graph_mem, bits_per_kmer;

  // reading stream: all covgs + edges
  // reading file: one bit per kmer per colour for 'in colour'
  bits_per_kmer = sizeof(BinaryKmer)*8;

  if(reading_stream) {
    bits_per_kmer += ncols * 8 * (sizeof(Edges) + sizeof(Covg));
  } else {
    bits_per_kmer += ncols; // in colour
  }

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
  int alloc_flags = reading_stream ? DBG_ALLOC_EDGES | DBG_ALLOC_COVGS
                                   : DBG_ALLOC_NODE_IN_COL;

  dBGraph db_graph;
  db_graph_alloc(&db_graph, file.hdr.kmer_size,
                 ncols, reading_stream ? ncols : 1,
                 kmers_in_hash, alloc_flags);

  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);

  // We need to load the graph for both --pop and --all since we need to check
  // if the next kmer is in each of the colours
  graph_load(&file, gprefs, NULL);

  if(add_pop_edges) status("Inferring edges from population...\n");
  else status("Inferring all missing edges...\n");

  size_t num_kmers_edited;

  if(reading_stream)
  {
    // Reading STDIN, writing STDOUT/file
    ctx_assert(fout != NULL);
    num_kmers_edited = infer_edges(num_of_threads, add_all_edges, &db_graph);
    graph_write_header(fout, &file.hdr);
    graph_write_all_kmers_direct(fout, &db_graph, false, &file.hdr);
  }
  else if(fout == NULL) {
    // Reading from file, writing to same file
    num_kmers_edited = inferedges_on_mmap(&db_graph, add_all_edges, &file);
  } else {
    // Reading from file, writing to STDOUT/file
    num_kmers_edited = inferedges_on_file(&db_graph, add_all_edges, &file, fout);
  }

  if(fout != NULL && fout != stdout) fclose(fout);

  char modified_str[100], kmers_str[100];
  ulong_to_str(num_kmers_edited, modified_str);
  ulong_to_str(hash_table_nkmers(&db_graph.ht), kmers_str);

  double modified_rate = 0;
  if(hash_table_nkmers(&db_graph.ht))
    modified_rate = (100.0 * num_kmers_edited) / hash_table_nkmers(&db_graph.ht);

  status("%s of %s (%.2f%%) nodes modified\n",
         modified_str, kmers_str, modified_rate);

  if(editing_file)
  {
    // Close and re-open
    fclose(file.fh);
    file.fh = NULL;
    futil_update_timestamp(file.fltr.path.b);
  }

  graph_file_close(&file);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
