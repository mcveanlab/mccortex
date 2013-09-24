#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "file_reader.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" inferedges [options] <pop.ctx>\n"
"  Infer edges in a population graph.  \n"
"\n"
"  -m <mem>   Memory to use (e.g. 100G or 12M)\n"
"  -k <kmer>  Kmer size\n";

// If two kmers are in a sample and the population has an edges between them,
// Add edge to sample

// Return 0 if no change
// Return 1 if changed
static inline int infer_edges(const BinaryKmer node_bkey, Edges *edges,
                              const Covg *covgs, dBGraph *db_graph)
{
  Edges uedges = 0, iedges = 0xf, add_edges, edge;
  size_t orient, nuc, col, kmer_size = db_graph->kmer_size;
  BinaryKmer bkey, bkmer;
  hkey_t next;

  Edges newedges[db_graph->num_of_cols];

  for(col = 0; col < db_graph->num_of_cols; col++) {
    uedges |= edges[col]; // union of edges
    iedges &= edges[col]; // intersection of edges
    newedges[col] = edges[col];
  }

  add_edges = uedges & ~iedges;

  if(!add_edges) return 0;

  for(orient = 0; orient < 2; orient++)
  {
    bkmer = node_bkey;
    if(orient == FORWARD) binary_kmer_left_shift_one_base(&bkmer, kmer_size);
    else binary_kmer_right_shift_one_base(&bkmer);

    for(nuc = 0; nuc < 4; nuc++)
    {
      edge = nuc_orient_to_edge(nuc, orient);
      if(add_edges & edge)
      {
        // get next bkmer, look up in graph
        if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
        else binary_kmer_set_first_nuc(&bkmer, binary_nuc_complement(nuc), kmer_size);

        bkey = db_node_get_key(bkmer, kmer_size);
        next = hash_table_find(&db_graph->ht, bkey);

        for(col = 0; col < db_graph->num_of_cols; col++)
          if(covgs[col] > 0 && db_node_has_col(db_graph, next, col))
            newedges[col] |= edge;
      }
    }
  }

  int cmp = memcmp(edges, newedges, sizeof(Edges)*db_graph->num_of_cols);
  memcpy(edges, newedges, sizeof(Edges)*db_graph->num_of_cols);

  return (cmp != 0);
}

int ctx_infer_edges(CmdArgs *args)
{
  cmd_accept_options(args, "mh");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 1) print_usage(usage, NULL);
  
  const char *path = argv[0];
  dBGraph db_graph;
  boolean is_binary = false;
  GraphFileHeader gheader = {.capacity = 0};

  if(!graph_file_probe(path, &is_binary, &gheader))
    print_usage(usage, "Cannot read input binary file: %s", path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", path);

  if(!test_file_writable(path))
    print_usage(usage, "Cannot write to file: %s", path);

  //
  // Decide on memory
  //
  size_t extra_bits_per_kmer = sizeof(Edges)*8 + gheader.num_of_cols;
  size_t kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                               gheader.num_of_kmers, true);

  db_graph_alloc(&db_graph, gheader.kmer_size,
                 gheader.num_of_cols, 1, kmers_in_hash);

  // In colour
  size_t words64_per_col = round_bits_to_words64(kmers_in_hash);
  db_graph.node_in_cols = calloc2(words64_per_col*gheader.num_of_cols, sizeof(uint64_t));
  db_graph.col_edges = calloc2(kmers_in_hash, sizeof(Edges));

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .db_graph = &db_graph,
                           .merge_colours = false,
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = false};

  graph_load(path, &prefs, stats, NULL);

  status("Inferring edges from population...\n");

  // Read again
  BinaryKmer bkmer;
  Edges edges[db_graph.num_of_cols];
  Covg covgs[db_graph.num_of_cols];

  FILE *fh = fopen(path, "r+");
  if(fh == NULL) die("Cannot open: %s", path);

  graph_file_read_header(fh, &gheader, true, path);

  size_t num_nodes_modified = 0;

  long edges_len = sizeof(Edges) * gheader.num_of_cols;
  while(graph_file_read_kmer(fh, &gheader, path, bkmer.b, covgs, edges))
  {
    if(infer_edges(bkmer, edges, covgs, &db_graph))
    {
      if(fseek(fh, -edges_len, SEEK_CUR) != 0) die("fseek error: %s", path);
      fwrite(edges, 1, edges_len, fh);
      num_nodes_modified++;
    }
  }

  fclose(fh);

  char modified_str[100], kmers_str[100];
  ulong_to_str(num_nodes_modified, modified_str);
  ulong_to_str(db_graph.ht.unique_kmers, kmers_str);
  status("%s of %s (%.2f%%) nodes modified\n", modified_str, kmers_str,
         (100.0 * num_nodes_modified) / db_graph.ht.unique_kmers);

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
