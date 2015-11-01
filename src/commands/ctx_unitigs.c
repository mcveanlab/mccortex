#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "supernode.h"
#include "graphs_load.h"
#include "gpath_checks.h"
#include "unitig_graph.h"

const char unitigs_usage[] =
"usage: "CMD" unitigs [options] <in.ctx> [<in2.ctx> ...]\n"
"\n"
"  Print unitigs with k-1 bases of overlap.\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -f, --force           Overwrite output files\n"
"  -o, --out <out.txt>   Save output graph file [default: STDOUT]\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
//
"  -F, --fasta           Print in FASTA format (default)\n"
"  -g, --gfa             Print in Graphical Fragment Assembly (GFA) format\n"
"  -d, --dot             Print in graphviz (DOT) format\n"
"  -P, --points          Used with --dot, print contigs as points\n"
"\n"
"  e.g. "CMD" unitigs --dot in.ctx | dot -Tpdf > in.pdf\n"
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
  {"fasta",        no_argument,       NULL, 'F'},
  {"gfa",          no_argument,       NULL, 'g'},
  {"dot",          no_argument,       NULL, 'd'},
  {"points",       no_argument,       NULL, 'P'},
  {NULL, 0, NULL, 0}
};


typedef enum {
  PRINT_FASTA = 0,
  PRINT_GFA = 1,
  PRINT_DOT = 2
} UnitigSyntax;

const char *syntax_strs[3] = {"FASTA", "GFA", "DOT (Graphviz)"};


typedef struct
{
  const dBGraph *db_graph;
  size_t nthreads;
  uint8_t *visited;
  UnitigSyntax syntax;
  FILE *fout;
  pthread_mutex_t outlock;
  UnitigKmerGraph ugraph;
  volatile size_t num_unitigs;
} UnitigPrinter;

/**
 * @param right_edge is true iff we this kmer is the last in a unitig
 */
static inline void _print_edge(hkey_t node, bool right_edge,
                               BinaryKmer bkey, Edges edges,
                               UnitigEnd uend0,
                               UnitigPrinter *p)
{
  // DOT: leave from east end if +, west end if -
  //      connect to west end if +, east end if -
  const char dot_exit[2] = "ew", dot_join[2] = "we", gfa_orient[2] = "+-";
  size_t i, n;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  Orientation orient = right_edge ? uend0.rorient : !uend0.lorient;
  // Unitig orientations
  Orientation ut_or0 = right_edge ? FORWARD : REVERSE, ut_or1;

  n = db_graph_next_nodes(p->db_graph, bkey, orient, edges,
                          next_nodes, next_nucs);

  for(i = 0; i < n; i++)
  {
    UnitigEnd uend1 = p->ugraph.unitig_ends[next_nodes[i].key];

    // Debugging
    if(!uend1.assigned) {
      char tmpstr[100];
      db_node_to_str(p->db_graph, next_nodes[i], tmpstr);
      status(" -> node %zu [%s]", (size_t)uend1.unitigid, tmpstr);
    }

    ctx_assert(next_nodes[i].key != HASH_NOT_FOUND);
    ctx_assert(uend1.assigned);

    ut_or1 = next_nodes[i].orient == uend1.lorient ? FORWARD : REVERSE;

    // Don't do reverse-to-reverse links when node links to itself,
    // these are duplicates of forward-to-forward
    if(node < next_nodes[i].key ||
       (node == next_nodes[i].key && ut_or0 + ut_or1 < 2))
    {
      pthread_mutex_lock(&p->outlock);

      switch(p->syntax) {
        case PRINT_DOT:
          fprintf(p->fout, "  node%zu:%c -> node%zu:%c\n",
                  (size_t)uend0.unitigid, dot_exit[ut_or0],
                  (size_t)uend1.unitigid, dot_join[ut_or1]);
          break;
        case PRINT_GFA:
          fprintf(p->fout, "L\tnode%zu\t%c\tnode%zu\t%c\t%zuM\n",
                  (size_t)uend0.unitigid, gfa_orient[ut_or0],
                  (size_t)uend1.unitigid, gfa_orient[ut_or1],
                  p->db_graph->kmer_size - 1);
          break;
        default: die("Bad syntax: %i", p->syntax);
      }

      pthread_mutex_unlock(&p->outlock);
    }
  }
}

// For every kmer in the graph, we run this function
static inline bool print_edges(hkey_t hkey, size_t threadid, void *arg)
{
  (void)threadid;
  UnitigPrinter *p = (UnitigPrinter*)arg;
  UnitigEnd uend0 = p->ugraph.unitig_ends[hkey];
  BinaryKmer bkey; Edges edges;

  // Check if node is an end of a supernode
  if(uend0.assigned) {
    bkey = db_node_get_bkmer(p->db_graph, hkey);
    edges = db_node_get_edges(p->db_graph, hkey, 0);

    if(uend0.left) {
      _print_edge(hkey, false, bkey, edges, uend0, p);
    }
    if(uend0.right) {
      _print_edge(hkey, true, bkey, edges, uend0, p);
    }
  }

  return false; // keep iterating
}

// Does not allocate ugraph
void unitig_printer_init(UnitigPrinter *printer, const dBGraph *db_graph,
                         size_t nthreads, UnitigSyntax syntax, FILE *fout)
{
  memset(printer, 0, sizeof(UnitigPrinter));
  printer->db_graph = db_graph;
  printer->syntax = syntax;
  printer->fout = fout;
  printer->nthreads = nthreads;
  printer->num_unitigs = 0;
  printer->visited = ctx_calloc(roundup_bits2bytes(db_graph->ht.capacity), 1);
  if(pthread_mutex_init(&printer->outlock, NULL) != 0) die("Mutex init failed");
}

void unitig_printer_destroy(UnitigPrinter *printer)
{
  pthread_mutex_destroy(&printer->outlock);
  unitig_graph_dealloc(&printer->ugraph);
  ctx_free(printer->visited);
}

static void print_unitig_fasta(dBNodeBuffer nbuf, size_t threadid, void *arg)
{
  (void)threadid;
  UnitigPrinter *p = (UnitigPrinter*)arg;
  pthread_mutex_lock(&p->outlock);
  size_t idx = p->num_unitigs++;
  fprintf(p->fout, ">supernode%zu\n", idx);
  db_nodes_print(nbuf.b, nbuf.len, p->db_graph, p->fout);
  fputc('\n', p->fout);
  pthread_mutex_unlock(&p->outlock);
}

static void print_unitig_dot(const dBNode *nodes, size_t num_nodes,
                             size_t unitig_idx, void *arg)
{
  UnitigPrinter *p = (UnitigPrinter*)arg;
  pthread_mutex_lock(&p->outlock);
  fprintf(p->fout, "  node%zu [label=", unitig_idx);
  db_nodes_print(nodes, num_nodes, p->db_graph, p->fout);
  fputs("]\n", p->fout);
  pthread_mutex_unlock(&p->outlock);
}

static void print_dot_syntax(UnitigPrinter *p, bool dot_use_points)
{
  fputs("digraph G {\n", p->fout);
  fputs("  edge [dir=both arrowhead=none arrowtail=none color=\"blue\"]\n", p->fout);
  fprintf(p->fout, "  node [%s, fontname=courier, fontsize=9]\n",
          dot_use_points ? "shape=point, label=none" : "shape=none");

  unitig_graph_create(&p->ugraph, p->nthreads, p->visited,
                      print_unitig_dot, p);

  p->num_unitigs = p->ugraph.num_unitigs;

  // Now print edges
  fputc('\n', p->fout);
  hash_table_iterate(&p->db_graph->ht, p->nthreads, print_edges, p);
  fputs("}\n", p->fout);
}


static void print_unitig_gfa(const dBNode *nodes, size_t num_nodes,
                             size_t unitig_idx, void *arg)
{
  UnitigPrinter *p = (UnitigPrinter*)arg;
  pthread_mutex_lock(&p->outlock);
  fprintf(p->fout, "S\tnode%zu\t", unitig_idx);
  db_nodes_print(nodes, num_nodes, p->db_graph, p->fout);
  fputc('\n', p->fout);
  pthread_mutex_unlock(&p->outlock);
}

static void print_gfa_syntax(UnitigPrinter *p)
{
  fputs("H VN:Z:1.0\n", p->fout);

  unitig_graph_create(&p->ugraph, p->nthreads, p->visited,
                      print_unitig_gfa, p);

  p->num_unitigs = p->ugraph.num_unitigs;
  // Now print edges
  hash_table_iterate(&p->db_graph->ht, p->nthreads, print_edges, p);
}

// Returns 0 on success, otherwise != 0
int ctx_unitigs(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  UnitigSyntax syntax = PRINT_FASTA;
  bool dot_use_points = false;

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
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 't': cmd_check(!nthreads, cmd); nthreads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'F': cmd_check(!syntax, cmd); syntax = PRINT_FASTA; break;
      case 'g': cmd_check(!syntax, cmd); syntax = PRINT_GFA; break;
      case 'd': cmd_check(!syntax, cmd); syntax = PRINT_DOT; break;
      case 'P': cmd_check(!dot_use_points, cmd); dot_use_points = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" unitigs -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(dot_use_points && syntax == PRINT_FASTA)
    cmd_print_usage("--point is only for use with --dot");

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;

  if(optind >= argc) cmd_print_usage(NULL);

  size_t i, num_gfiles = (size_t)(argc - optind);
  char **gfile_paths = argv + optind;

  if(dot_use_points && syntax != PRINT_DOT)
    cmd_print_usage("--points only valid with --graphviz / --dot");

  ctx_assert(num_gfiles > 0);

  // Open graph files
  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

  graph_files_open(gfile_paths, gfiles, num_gfiles,
                   &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 + 1;
  if(syntax != PRINT_FASTA) bits_per_kmer += sizeof(UnitigEnd) * 8;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  status("Output in %s format to %s\n", syntax_strs[syntax],
         futil_outpath_str(out_path));

  //
  // Open output file
  //

  // Print to stdout unless --out <out> is specified
  FILE *fout = futil_fopen_create(out_path, "w");

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, 1, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES);

  UnitigPrinter printer;
  unitig_printer_init(&printer, &db_graph, nthreads, syntax, fout);

  if(syntax == PRINT_DOT || syntax == PRINT_GFA)
    unitig_graph_alloc(&printer.ugraph, &db_graph);

  // Load graphs
  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);

  for(i = 0; i < num_gfiles; i++) {
    file_filter_flatten(&gfiles[i].fltr, 0);
    graph_load(&gfiles[i], gprefs, NULL);
    graph_file_close(&gfiles[i]);
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  switch(syntax)
  {
    case PRINT_FASTA:
      status("Printing unitgs in FASTA using %zu threads", nthreads);
      supernodes_iterate(nthreads, printer.visited, &db_graph,
                         print_unitig_fasta, &printer);
      break;
    case PRINT_GFA:
      print_gfa_syntax(&printer);
      break;
    case PRINT_DOT:
      print_dot_syntax(&printer, dot_use_points);
      break;
    default:
      die("Invalid print syntax: %i", syntax);
  }

  char num_unitigs_str[50];
  ulong_to_str(printer.num_unitigs, num_unitigs_str);
  status("Dumped %s unitigs\n", num_unitigs_str);

  fclose(fout);

  unitig_printer_destroy(&printer);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
