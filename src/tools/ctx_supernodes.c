#include "global.h"

#include "tools.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "binary_kmer.h"
#include "supernode.h"

const char supernodes_usage[] =
"usage: "CMD" supernodes [options] <in.ctx> [<in2.ctx> ...]\n"
"  Print supernodes with k-1 bases of overlap.\n"
"\n"
"  Options:  --memory <M> Memory to use (e.g. 2GB)\n"
"            --out <out>  specify output [default: STDOUT]\n"
"            --graphviz   print in graphviz format\n"
"            --point      with --graphviz, print contigs as points\n"
"\n"
"  e.g. ctx31 supernodes --graphviz in.ctx | dot -Tpdf > in.pdf\n";

typedef struct {
  size_t nodeid:57, left:1, right:1, lorient:1, rorient:1, assigned:1;
} sndata_t;


#define PRINT_FASTA 0
#define PRINT_DOT 1

static size_t supernode_idx = 0;

// Store ends of supernode currently stored in `nodes` and `orients` arrays
static inline void dot_store_ends(const dBNodeBuffer *nbuf, sndata_t *supernodes)
{
  ctx_assert(supernodes[nbuf->data[0].key].assigned == 0);
  ctx_assert(supernodes[nbuf->data[nbuf->len-1].key].assigned == 0);

  sndata_t supernode0 = {.nodeid = supernode_idx, .assigned = 1,
                         .left = 1, .right = (nbuf->len == 1),
                         .lorient = nbuf->data[0].orient,
                         .rorient = nbuf->data[nbuf->len-1].orient};

  sndata_t supernode1 = {.nodeid = supernode_idx, .assigned = 1,
                         .left = (nbuf->len == 1), .right = 1,
                         .lorient = nbuf->data[0].orient,
                         .rorient = nbuf->data[nbuf->len-1].orient};

  supernodes[nbuf->data[0].key] = supernode0;
  supernodes[nbuf->data[nbuf->len-1].key] = supernode1;
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

  for(i = 0; i < n; i++)
  {
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
static inline void dot_print_edges(hkey_t node, sndata_t *supernodes,
                                   FILE *fout, const dBGraph *db_graph)
{
  sndata_t snode0 = supernodes[node];
  BinaryKmer bkmer; Edges edges;

  // Check if node is an end of a supernode
  if(snode0.assigned) {
    bkmer = db_node_get_bkmer(db_graph, node);
    edges = db_node_get_edges(db_graph, 0, node);

    if(snode0.left) {
      dot_print_edges2(node, bkmer, edges, !snode0.lorient, snode0,
                       supernodes, fout, db_graph);
    }
    if(snode0.right) {
      dot_print_edges2(node, bkmer, edges, snode0.rorient, snode0,
                       supernodes, fout, db_graph);
    }
  }
}

static void dump_supernodes(hkey_t hkey, FILE *fout, int print_syntax,
                            dBNodeBuffer *nbuf, sndata_t *supernodes,
                            uint64_t *visited, const dBGraph *db_graph)
{
  ctx_assert(print_syntax == PRINT_FASTA || supernodes != NULL);

  size_t i;

  if(!bitset_get(visited, hkey))
  {
    db_node_buf_reset(nbuf);
    supernode_find(hkey, nbuf, db_graph);
    for(i = 0; i < nbuf->len; i++) bitset_set(visited, nbuf->data[i].key);

    supernode_normalise(nbuf->data, nbuf->len, db_graph);

    switch(print_syntax) {
      case PRINT_FASTA:
        fprintf(fout, ">supernode%zu\n", supernode_idx);
        db_nodes_print(nbuf->data, nbuf->len, db_graph, stdout);
        fputc('\n', fout);
        break;
      case PRINT_DOT:
        dot_store_ends(nbuf, supernodes);
        fprintf(fout, "  node%zu [label=", supernode_idx);
        db_nodes_print(nbuf->data, nbuf->len, db_graph, stdout);
        fputs("]\n", fout);
        break;
    }

    supernode_idx++;
  }
}

static void dump_dot_syntax(FILE *fout, int print_syntax, boolean dot_use_points,
                            dBNodeBuffer *nbuf, uint64_t *visited,
                            dBGraph *db_graph)
{
  fputs("digraph G {\n", fout);
  fputs("  edge [dir=both arrowhead=none arrowtail=none color=\"blue\"]\n", fout);
  fprintf(fout, "  node [%s, fontname=courier, fontsize=9]\n",
          dot_use_points ? "shape=point, label=none" : "shape=none");

  sndata_t *supernodes = calloc2(db_graph->ht.capacity, sizeof(sndata_t));

  HASH_ITERATE(&db_graph->ht, dump_supernodes,
               fout, print_syntax, nbuf, supernodes, visited, db_graph);

  // Now print edges
  fprintf(fout, "\n");
  HASH_ITERATE(&db_graph->ht, dot_print_edges, supernodes, fout, db_graph);

  free(supernodes);
  fputs("}\n", fout);
}

// Returns 0 on success, otherwise != 0
int ctx_supernodes(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;

  size_t i, num_files;
  char **paths;
  uint64_t max_ctx_kmers = 0;
  boolean dot_use_points = false;
  int print_syntax = PRINT_FASTA;

  while(argc > 0 && argv[0][0] == '-') {
    if(!strcasecmp(argv[0],"--dot") || !strcasecmp(argv[0],"--graphviz")) {
      print_syntax = PRINT_DOT; argv++; argc--;
    }
    else if(!strcasecmp(argv[0],"--points")) {
      dot_use_points = true; argv++; argc--;
    }
    else cmd_print_usage("Unknown argument: %s", argv[0]);
  }

  num_files = (size_t)argc;
  paths = argv;

  if(num_files == 0) cmd_print_usage(NULL);

  if(dot_use_points && print_syntax != PRINT_DOT)
    cmd_print_usage("--points only valid with --graphviz / --dot");

  GraphFileReader files[num_files];

  for(i = 0; i < num_files; i++)
  {
    files[i] = INIT_GRAPH_READER;
    graph_file_open(&files[i], paths[i], true);

    if(files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, files[i].hdr.kmer_size);
    }

    max_ctx_kmers = MAX2(max_ctx_kmers, files[i].hdr.num_of_kmers);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  bits_per_kmer = (sizeof(Edges) + sizeof(sndata_t)*(print_syntax==PRINT_DOT))*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, max_ctx_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  const char *syntax_str[3] = {"FASTA", "DOT (Graphviz)"};
  status("Output in %s format to %s\n", syntax_str[print_syntax],
         args->output_file_set ? args->output_file : "STDOUT");

  //
  // Open output file
  //

  // Print to stdout unless --out <out> is specified
  FILE *fout = stdout;

  if(args->output_file_set) {
    if(strcmp(args->output_file, "-") == 0)
      args->output_file_set = false;
    else {
      fout = fopen(args->output_file, "w");
      if(fout == NULL) die("Cannot open output file: %s", args->output_file);
      status("Writing supernodes to %s\n", args->output_file);
    }
  }

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));

  // Visited
  size_t numwords64 = roundup_bits2words64(db_graph.ht.capacity);
  uint64_t *visited = calloc2(numwords64, sizeof(uint64_t));

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = true,
                              .must_exist_in_graph = false,
                              .empty_colours = false};

  for(i = 0; i < num_files; i++) {
    files[i].fltr.flatten = true;
    // files[i].fltr.intocol = 0;
    file_filter_update_intocol(&files[i].fltr, 0);
    graph_load(&files[i], gprefs, NULL);
  }

  hash_table_print_stats(&db_graph.ht);

  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 2048);

  // dump supernodes
  switch(print_syntax) {
    case PRINT_DOT:
      dump_dot_syntax(fout, print_syntax, dot_use_points,
                      &nbuf, visited, &db_graph);
      break;
    case PRINT_FASTA:
      HASH_ITERATE(&db_graph.ht, dump_supernodes,
                   fout, print_syntax, &nbuf, NULL, visited, &db_graph);
      break;
  }

  status("Dumped %zu supernodes\n", supernode_idx);

  if(args->output_file_set) fclose(fout);

  db_node_buf_dealloc(&nbuf);
  free(visited);
  free(db_graph.col_edges);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);

  return EXIT_SUCCESS;
}
