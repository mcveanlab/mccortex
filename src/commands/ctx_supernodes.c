#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "binary_kmer.h"
#include "supernode.h"

const char supernodes_usage[] =
"usage: "CMD" supernodes [options] <in.ctx> [<in2.ctx> ...]\n"
"\n"
"  Print supernodes with k-1 bases of overlap.\n"
"\n"
"  -h, --help            This help message\n"
"  -f, --force           Overwrite output files\n"
"  -o, --out <out.txt>   Save output graph file [default: STDOUT]\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>       Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
//
"  -d, --dot             Print in graphviz (DOT) format\n"
"  -p, --points          Used with --dot, print contigs as points\n"
"\n"
"  e.g. ctx31 supernodes --dot in.ctx | dot -Tpdf > in.pdf\n"
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
  {"graphviz",     no_argument,       NULL, 'g'}, // obsolete: use dot
  {"dot",          no_argument,       NULL, 'd'},
  {"points",       no_argument,       NULL, 'p'},
  {NULL, 0, NULL, 0}
};

// Each supernode is packed into 64 bits, and each kmer has supernode info
typedef struct {
  size_t nodeid:57, left:1, right:1, lorient:1, rorient:1, assigned:1;
} sndata_t;

#define PRINT_FASTA 0
#define PRINT_DOT 1

// Store ends of supernode currently stored in `nodes` and `orients` arrays
static inline void dot_store_ends(size_t snidx, dBNodeBuffer nbuf,
                                  sndata_t *supernodes)
{
  ctx_assert(supernodes[nbuf.data[0].key].assigned == 0);
  ctx_assert(supernodes[nbuf.data[nbuf.len-1].key].assigned == 0);

  sndata_t supernode0 = {.nodeid = snidx, .assigned = 1,
                         .left = 1, .right = (nbuf.len == 1),
                         .lorient = nbuf.data[0].orient,
                         .rorient = nbuf.data[nbuf.len-1].orient};

  sndata_t supernode1 = {.nodeid = snidx, .assigned = 1,
                         .left = (nbuf.len == 1), .right = 1,
                         .lorient = nbuf.data[0].orient,
                         .rorient = nbuf.data[nbuf.len-1].orient};

  supernodes[nbuf.data[0].key] = supernode0;
  supernodes[nbuf.data[nbuf.len-1].key] = supernode1;
}

static inline void dot_print_edges2(hkey_t node, BinaryKmer bkmer, Edges edges,
                                    Orientation orient, sndata_t snode0,
                                    sndata_t *supernodes,
                                    FILE *fout,
                                    const dBGraph *db_graph)
{
  size_t i, n, side0, side1;
  dBNode next_nodes[4];
  Nucleotide next_nucs[4];
  sndata_t snode1;
  const char coords[2] = "we";

  n = db_graph_next_nodes(db_graph, bkmer, orient, edges,
                          next_nodes, next_nucs);

  // side0 = snode0.left && snode0.right ? !orient : snode0.right;
  // fprintf(fout, "node%zu:%c\n", (size_t)snode0.nodeid, coords[side0]);

  // char tmp[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(bkmer, db_graph->kmer_size, tmp);
  // printf("bkmer: %s [orient=%i, n=%zu]\n", tmp, orient, n);

  for(i = 0; i < n; i++)
  {
    // printf("  i=%zu\n", i);
    snode1 = supernodes[next_nodes[i].key];
    ctx_assert(next_nodes[i].key != HASH_NOT_FOUND);
    ctx_assert(snode1.assigned);

    // if left && right then supernode is 1 kmer long
    // orient: FORWARD == 0, 'e'; REVERSE == 1, 'w'; => !orient
    // next_orients: need to negate twice so just next_nodes[i].orient
    side0 = snode0.left && snode0.right ? !orient : snode0.right;
    side1 = snode1.left && snode1.right ? next_nodes[i].orient : snode1.right;

    if(node < next_nodes[i].key || (node == next_nodes[i].key && side0 <= side1)) {
      fprintf(fout, "  node%zu:%c -> node%zu:%c\n",
              (size_t)snode0.nodeid, coords[side0],
              (size_t)snode1.nodeid, coords[side1]);
    }
  }
}

// For every kmer in the graph, we run this function
static inline void dot_print_edges(hkey_t hkey, sndata_t *supernodes,
                                   FILE *fout, const dBGraph *db_graph)
{
  sndata_t snode0 = supernodes[hkey];
  BinaryKmer bkmer; Edges edges;

  // Check if node is an end of a supernode
  if(snode0.assigned) {
    bkmer = db_node_get_bkmer(db_graph, hkey);
    edges = db_node_get_edges(db_graph, hkey, 0);

    if(snode0.left) {
      dot_print_edges2(hkey, bkmer, edges, !snode0.lorient, snode0,
                       supernodes, fout, db_graph);
    }
    if(snode0.right) {
      dot_print_edges2(hkey, bkmer, edges, snode0.rorient, snode0,
                       supernodes, fout, db_graph);
    }
  }
}

struct SupernodePrinter
{
  FILE *fout;
  pthread_mutex_t outlock;
  sndata_t *supernodes;
  size_t *supernode_idx;
  const int print_syntax;
  const dBGraph *db_graph;
};

static void print_supernodes(dBNodeBuffer nbuf, size_t threadid, void *arg)
{
  (void)threadid;
  struct SupernodePrinter *prtr = (struct SupernodePrinter*)arg;

  supernode_normalise(nbuf.data, nbuf.len, prtr->db_graph);

  size_t idx = __sync_fetch_and_add((size_t volatile*)prtr->supernode_idx, 1);
  FILE *fout = prtr->fout;

  if(prtr->print_syntax == PRINT_DOT)
    dot_store_ends(idx, nbuf, prtr->supernodes);

  pthread_mutex_lock(&prtr->outlock);

  if(prtr->print_syntax == PRINT_FASTA) {
    fprintf(fout, ">supernode%zu\n", idx);
    db_nodes_print(nbuf.data, nbuf.len, prtr->db_graph, fout);
    fputc('\n', fout);
  }
  else {
    ctx_assert(prtr->print_syntax == PRINT_DOT);
    fprintf(fout, "  node%zu [label=", idx);
    db_nodes_print(nbuf.data, nbuf.len, prtr->db_graph, fout);
    fputs("]\n", fout);
  }

  pthread_mutex_unlock(&prtr->outlock);
}

// Returns number of supernodes printed
static size_t print_all_supernodes(size_t nthreads, FILE *fout,
                                   int print_syntax, sndata_t *supernodes,
                                   uint8_t *visited, const dBGraph *db_graph)
{
  ctx_assert(print_syntax == PRINT_FASTA || supernodes != NULL);

  size_t next_snode_idx = 0;
  struct SupernodePrinter printer = {.fout = fout,
                                     .supernodes = supernodes,
                                     .supernode_idx = &next_snode_idx,
                                     .print_syntax = print_syntax,
                                     .db_graph = db_graph};

  if(pthread_mutex_init(&printer.outlock, NULL) != 0) die("Mutex init failed");
  supernodes_iterate(nthreads, visited, db_graph, print_supernodes, &printer);
  pthread_mutex_destroy(&printer.outlock);

  return next_snode_idx;
}

static size_t print_dot_syntax(size_t nthreads, FILE *fout,
                               int print_syntax, bool dot_use_points,
                               uint8_t *visited, const dBGraph *db_graph)
{
  fputs("digraph G {\n", fout);
  fputs("  edge [dir=both arrowhead=none arrowtail=none color=\"blue\"]\n", fout);
  fprintf(fout, "  node [%s, fontname=courier, fontsize=9]\n",
          dot_use_points ? "shape=point, label=none" : "shape=none");

  sndata_t *supernodes = ctx_calloc(db_graph->ht.capacity, sizeof(sndata_t));

  size_t num_snodes = print_all_supernodes(nthreads, fout, print_syntax,
                                           supernodes, visited, db_graph);

  // Now print edges
  fputc('\n', fout);
  HASH_ITERATE(&db_graph->ht, dot_print_edges, supernodes, fout, db_graph);
  fputs("}\n", fout);

  ctx_free(supernodes);

  return num_snodes;
}

// Returns 0 on success, otherwise != 0
int ctx_supernodes(int argc, char **argv)
{
  size_t nthreads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  int print_syntax = PRINT_FASTA;
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
      case 'g': // --graphviz is the same as --dot, drop through case
      case 'd': if(print_syntax) die("%s set twice", cmd); print_syntax=PRINT_DOT; break;
      case 'p': if(dot_use_points) die("%s set twice", cmd); dot_use_points=true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" supernodes -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(dot_use_points && print_syntax == PRINT_FASTA)
    cmd_print_usage("--point is only for use with --dot");

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(nthreads == 0) nthreads = DEFAULT_NTHREADS;

  if(optind >= argc) cmd_print_usage(NULL);

  size_t i, num_gfiles = (size_t)(argc - optind);
  char **gfile_paths = argv+optind;

  if(dot_use_points && print_syntax != PRINT_DOT)
    cmd_print_usage("--points only valid with --graphviz / --dot");

  ctx_assert(num_gfiles > 0);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

  graph_files_open(gfile_paths, gfiles, num_gfiles,
                   &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 + 1;
  if(print_syntax == PRINT_DOT) bits_per_kmer += sizeof(sndata_t) * 8;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  const char *syntax_str[3] = {"FASTA", "DOT (Graphviz)"};
  status("Output in %s format to %s\n", syntax_str[print_syntax],
         futil_outpath_str(out_path));

  //
  // Open output file
  //

  // Print to stdout unless --out <out> is specified
  FILE *fout = futil_open_output(out_path);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));

  uint8_t *visited = ctx_calloc(roundup_bits2bytes(db_graph.ht.capacity), 1);

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = false};

  for(i = 0; i < num_gfiles; i++) {
    gfiles[i].fltr.flatten = true;
    file_filter_update_intocol(&gfiles[i].fltr, 0);
    graph_load(&gfiles[i], gprefs, NULL);
    graph_file_close(&gfiles[i]);
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  status("Printing supernodes using %zu threads", nthreads);
  size_t num_snodes;

  // dump supernodes
  if(print_syntax == PRINT_FASTA) {
    num_snodes = print_all_supernodes(nthreads, fout,
                                      print_syntax, NULL,
                                      visited, &db_graph);
  }
  else {
    num_snodes = print_dot_syntax(nthreads, fout,
                                  print_syntax, dot_use_points,
                                  visited, &db_graph);
  }

  char num_snodes_str[50];
  ulong_to_str(num_snodes, num_snodes_str);
  status("Dumped %s supernodes\n", num_snodes_str);

  fclose(fout);

  ctx_free(visited);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
