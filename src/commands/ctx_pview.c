#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "gpath_reader.h"
#include "gpath_checks.h"

const char pview_usage[] =
"usage: "CMD" pview [options] [-p <in.ctp>] <in.ctx> [in2.ctx ...]\n"
"\n"
"  View cortex path files (.ctp).\n"
"\n"
"  -h, --help             This help message\n"
"  -q, --quiet            Silence status output normally printed to STDERR\n"
"  -m, --memory <mem>     Memory to use\n"
"  -n, --nkmers <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -p, --paths <in.ctp>   Load path file (can specify multiple times)\n"
"\n";


static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"paths",        required_argument, NULL, 'p'},
  {NULL, 0, NULL, 0}
};

static int _print_paths(hkey_t hkey, dBNodeBuffer *nbuf, FILE *fout,
                        const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;
  const GPath *gpath, *first_gpath = gpath_store_fetch_traverse(gpstore, hkey);
  dBNode node = {.key = hkey, .orient = FORWARD};
  size_t i, npaths = 0, col, ncols = gpstore->gpset.ncols;

  for(gpath = first_gpath; gpath != NULL; gpath = gpath->next) npaths++;

  if(npaths == 0) return 0;

  BinaryKmer bkey = db_node_get_bkmer(db_graph, node.key);
  char bkeystr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(bkey, db_graph->kmer_size, bkeystr);
  fprintf(fout, "%s %zu\n", bkeystr, npaths);

  char dir[2];
  dir[FORWARD] = 'F';
  dir[REVERSE] = 'R';

  for(gpath = first_gpath; gpath != NULL; gpath = gpath->next)
  {
    // Print gpath
    node.orient = gpath->orient;

    size_t klen = gpath_set_get_klen(gpset, gpath);
    size_t njuncs = gpath->num_juncs;
    const uint8_t *nseen = gpath_set_get_nseen(gpset, gpath);

    // Find a colour this path is in
    for(col = 0; col < ncols && !gpath_has_colour(gpath, ncols, col); col++) {}
    if(col == ncols) die("path is not in any colours");

    db_node_buf_reset(nbuf);
    gpath_fetch(node, gpath, nbuf, col, db_graph);
    ctx_assert(nbuf->len == klen);

    fprintf(fout, "%c %zu %zu %zu", dir[gpath->orient], klen, njuncs,
                                    (size_t)nseen[0]);

    for(i = 1; i < ncols; i++)
      fprintf(fout, ",%zu", (size_t)nseen[i]);

    fputc(' ', fout);
    binary_seq_print(gpath->seq, njuncs, fout);
    fputc(' ', fout);
    db_nodes_print(nbuf->data, nbuf->len, db_graph, fout);
    fputc('\n', fout);
  }

  return 0; // 0 => keep iterating
}

int ctx_pview(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);


  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_add(&gpfiles, tmp_gpfile);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" clean -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(optind >= argc)   cmd_print_usage("Please give input graph files");
  if(gpfiles.len == 0) cmd_print_usage("Please give input path files");

  // Use remaining args as graph files
  char **gfile_paths = argv + optind;
  size_t i, num_gfiles = (size_t)(argc - optind);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(gfile_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, gpfiles.data, gpfiles.len, -1);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 +
                  sizeof(GPath*)*8 + ncols;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.data, gpfiles.len, ncols, rem_mem, true);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.data, gpfiles.len, path_mem, true, &db_graph);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_load(&gpfiles.data[i], GPATH_DIE_MISSING_KMERS, &db_graph);

  // Print paths
  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 1024);
  HASH_ITERATE(&db_graph.ht, _print_paths, &nbuf, stdout, &db_graph);
  db_node_buf_dealloc(&nbuf);

  // Close input path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_close(&gpfiles.data[i]);
  gpfile_buf_dealloc(&gpfiles);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
