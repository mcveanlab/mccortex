#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "graph_file_filter.h"

static const char usage[] =
"usage: "CMD" healthcheck <out.ctx>\n"
"  Load a graph into memory to check it is valid.\n"
"\n"
"  Options:\n"
"   -m <mem>      Memory to use\n"
"   -h <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n";

static inline void check_node(hkey_t node, const dBGraph *db_graph)
{
  Edges edges = db_node_col_edges_union(db_graph, node);
  BinaryKmer bkmer = db_node_bkmer(db_graph, node);
  size_t nfw_edges, nrv_edges;
  hkey_t nodes[8];
  Orientation orients[8];
  Nucleotide nucs[8];

  nfw_edges = db_graph_next_nodes(db_graph, bkmer, FORWARD, edges,
                                  nodes, orients, nucs);
  nrv_edges = db_graph_next_nodes(db_graph, bkmer, REVERSE, edges,
                                  nodes, orients, nucs);

  if((unsigned)__builtin_popcount(edges) < nfw_edges + nrv_edges) {
    char seq[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph->kmer_size, seq);
    die("Excess edges on node: %s [%zu,%zu]", seq, nfw_edges, nrv_edges);
  }
}

int ctx_health_check(CmdArgs *args)
{
  cmd_accept_options(args, "mh");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 1) print_usage(usage, NULL);

  char *ctx_path = argv[0];
  GraphFileReader file = INIT_GRAPH_READER;
  graph_file_open(&file, ctx_path, true); // true => errors are fatal
  size_t ncols = graph_file_outncols(&file);

  // Decide on memory
  size_t extra_bits_per_kmer, kmers_in_hash;
  extra_bits_per_kmer = sizeof(Edges) * ncols * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        file.hdr.num_of_kmers, true);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, file.hdr.kmer_size, ncols, ncols, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity * ncols, sizeof(Edges));

  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = true};

  graph_load(&file, &prefs, NULL);

  HASH_TRAVERSE(&db_graph.ht, check_node, &db_graph);

  status("All looks good!");

  free(db_graph.col_edges);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
