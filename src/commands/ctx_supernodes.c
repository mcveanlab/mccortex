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
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -o, --out <out.ctx>   Save output graph file [default: STDOUT]\n"
"  -d, --graphviz        Print in graphviz (DOT) format\n"
"  -i, --point      with --graphviz, print contigs as points\n"
"\n"
"  e.g. ctx31 supernodes --graphviz in.ctx | dot -Tpdf > in.pdf\n"
"\n";

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
        db_nodes_print(nbuf->data, nbuf->len, db_graph, fout);
        fputc('\n', fout);
        break;
      case PRINT_DOT:
        dot_store_ends(nbuf, supernodes);
        fprintf(fout, "  node%zu [label=", supernode_idx);
        db_nodes_print(nbuf->data, nbuf->len, db_graph, fout);
        fputs("]\n", fout);
        break;
    }

    supernode_idx++;
  }
}

static void dump_dot_syntax(FILE *fout, int print_syntax, bool dot_use_points,
                            dBNodeBuffer *nbuf, uint64_t *visited,
                            dBGraph *db_graph)
{
  fputs("digraph G {\n", fout);
  fputs("  edge [dir=both arrowhead=none arrowtail=none color=\"blue\"]\n", fout);
  fprintf(fout, "  node [%s, fontname=courier, fontsize=9]\n",
          dot_use_points ? "shape=point, label=none" : "shape=none");

  sndata_t *supernodes = ctx_calloc(db_graph->ht.capacity, sizeof(sndata_t));

  HASH_ITERATE(&db_graph->ht, dump_supernodes,
               fout, print_syntax, nbuf, supernodes, visited, db_graph);

  // Now print edges
  fprintf(fout, "\n");
  HASH_ITERATE(&db_graph->ht, dot_print_edges, supernodes, fout, db_graph);

  ctx_free(supernodes);
  fputs("}\n", fout);
}

// Returns 0 on success, otherwise != 0
int ctx_supernodes(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;

  size_t i;
  char **paths;
  bool dot_use_points = false;
  int print_syntax = PRINT_FASTA;

  while(argc > 0 && argv[0][0] == '-' && argv[0][1]) {
    if(!strcmp(argv[0],"--dot") || !strcmp(argv[0],"-d") ||
       !strcmp(argv[0],"--graphviz"))
    {
      print_syntax = PRINT_DOT; argv++; argc--;
    }
    else if(!strcmp(argv[0],"--points") || !strcmp(argv[0],"--point") ||
            !strcmp(argv[0],"-i"))
    {
      dot_use_points = true; argv++; argc--;
    }
    else cmd_print_usage("Unknown argument: %s", argv[0]);
  }

  if(argc == 0) cmd_print_usage(NULL);

  const size_t num_gfiles = (size_t)argc;
  paths = argv;

  if(dot_use_points && print_syntax != PRINT_DOT)
    cmd_print_usage("--points only valid with --graphviz / --dot");

  ctx_assert(num_gfiles > 0);

  GraphFileReader gfiles[num_gfiles];
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;

  graph_files_open(paths, gfiles, num_gfiles,
                   &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  bits_per_kmer = (sizeof(Edges) + sizeof(sndata_t)*(print_syntax==PRINT_DOT))*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(args->mem_to_use, graph_mem);

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
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));

  // Visited
  size_t numwords64 = roundup_bits2words64(db_graph.ht.capacity);
  uint64_t *visited = ctx_calloc(numwords64, sizeof(uint64_t));

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

  hash_table_print_stats(&db_graph.ht);

  status("Printing supernodes...");

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
  ctx_free(visited);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
