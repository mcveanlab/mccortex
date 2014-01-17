#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "loading_stats.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" inferedges [options] <pop.ctx>\n"
"  Infer edges in a population graph.  \n"
"\n"
"  -m <mem>   Memory to use (e.g. 100G or 12M)\n"
"  -k <kmer>  Kmer size\n"
"  --pop      Add edges that are in the union only\n"
"  --all      Add all edges [default]\n";

// If two kmers are in a sample and the population has an edges between them,
// Add edge to sample

// Return 1 if changed; 0 otherwise
static inline int infer_pop_edges(const BinaryKmer node_bkey, Edges *edges,
                                  const Covg *covgs, dBGraph *db_graph)
{
  Edges uedges = 0, iedges = 0xf, add_edges, edge;
  size_t orient, nuc, col, kmer_size = db_graph->kmer_size;
  const size_t ncols = db_graph->num_of_cols;
  BinaryKmer bkey, bkmer;
  hkey_t next;
  Edges newedges[ncols];

  // char tmp[MAX_KMER_SIZE+1];
  // binary_kmer_to_str(node_bkey, db_graph->kmer_size, tmp);
  // status("Inferring %s", tmp);

  for(col = 0; col < ncols; col++) {
    uedges |= edges[col]; // union of edges
    iedges &= edges[col]; // intersection of edges
    newedges[col] = edges[col];
  }

  add_edges = uedges & ~iedges;

  if(!add_edges) return 0;

  for(orient = 0; orient < 2; orient++)
  {
    bkmer = (orient == FORWARD ? binary_kmer_left_shift_one_base(node_bkey, kmer_size)
                               : binary_kmer_right_shift_one_base(node_bkey));

    for(nuc = 0; nuc < 4; nuc++)
    {
      edge = nuc_orient_to_edge(nuc, orient);
      if(add_edges & edge)
      {
        // get next bkmer, look up in graph
        if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
        else binary_kmer_set_first_nuc(&bkmer, dna_nuc_complement(nuc), kmer_size);

        bkey = db_node_get_key(bkmer, kmer_size);
        next = hash_table_find(&db_graph->ht, bkey);
        assert(next != HASH_NOT_FOUND);

        for(col = 0; col < ncols; col++)
          if(covgs[col] > 0 && db_node_has_col(db_graph, next, col))
            newedges[col] |= edge;
      }
    }
  }

  int cmp = memcmp(edges, newedges, sizeof(Edges)*ncols);
  memcpy(edges, newedges, sizeof(Edges)*ncols);
  return (cmp != 0);
}

// Return 1 if changed; 0 otherwise
static inline int infer_all_edges(const BinaryKmer node_bkey, Edges *edges,
                                  const Covg *covgs, dBGraph *db_graph)
{
  Edges iedges = 0xff, edge;
  size_t orient, nuc, col, kmer_size = db_graph->kmer_size;
  const size_t ncols = db_graph->num_of_cols;
  BinaryKmer bkey, bkmer;
  hkey_t next;

  Edges newedges[ncols];
  memcpy(newedges, edges, ncols * sizeof(Edges));

  // intersection of edges
  for(col = 0; col < ncols; col++) iedges &= edges[col];

  for(orient = 0; orient < 2; orient++)
  {
    bkmer = (orient == FORWARD ? binary_kmer_left_shift_one_base(node_bkey, kmer_size)
                               : binary_kmer_right_shift_one_base(node_bkey));

    for(nuc = 0; nuc < 4; nuc++)
    {
      edge = nuc_orient_to_edge(nuc, orient);
      if(!(iedges & edge))
      {
        // edges are missing from some samples
        if(orient == FORWARD) binary_kmer_set_last_nuc(&bkmer, nuc);
        else binary_kmer_set_first_nuc(&bkmer, dna_nuc_complement(nuc), kmer_size);

        bkey = db_node_get_key(bkmer, kmer_size);
        next = hash_table_find(&db_graph->ht, bkey);

        if(next != HASH_NOT_FOUND) {
          for(col = 0; col < ncols; col++) {
            if(covgs[col] > 0 && db_node_has_col(db_graph, next, col)) {
              newedges[col] |= edge;
            }
          }
        }
      }
    }
  }

  int cmp = memcmp(edges, newedges, sizeof(Edges)*ncols);
  memcpy(edges, newedges, sizeof(Edges)*ncols);
  return (cmp != 0);
}


int ctx_infer_edges(CmdArgs *args)
{
  cmd_accept_options(args, "mn", usage);
  int argc = args->argc;
  char **argv = args->argv;
  // if(argc != 1) print_usage(usage, NULL);

  boolean add_pop_edges = false, add_all_edges = false;
  while(argc > 0 && argv[0][0] == '-') {
    if(strcmp(argv[0],"--all") == 0) {
      if(add_pop_edges) print_usage(usage, "Specify only one of --all --pop");
      argc--; argv++; add_all_edges = true;
    }
    else if(strcmp(argv[0],"--pop") == 0) {
      if(add_all_edges) print_usage(usage, "Specify only one of --all --pop");
      argc--; argv++; add_pop_edges = true;
    }
  }

  // Default to adding all edges
  if(!add_pop_edges && !add_all_edges) add_all_edges = true;

  if(argc != 1) print_usage(usage, NULL);

  char *path = argv[0];
  dBGraph db_graph;
  GraphFileReader file = INIT_GRAPH_READER;
  graph_file_open2(&file, path, true, "r+");

  if(!file.fltr.nofilter)
    print_usage(usage, "Inferedges with filter not implemented - sorry");

  if(!futil_is_file_writable(path))
    print_usage(usage, "Cannot write to file: %s", path);

  assert(file.hdr.num_of_cols == file.fltr.ncols);

  // file.fltr.intocol = 0;
  file_filter_update_intocol(&file.fltr, 0);
  file.fltr.flatten = false;

  //
  // Decide on memory
  //
  size_t kmers_in_hash, extra_bits_per_kmer = file.hdr.num_of_cols, graph_mem;
  kmers_in_hash = cmd_get_kmers_in_hash(args, extra_bits_per_kmer,
                                        file.hdr.num_of_kmers, true, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  db_graph_alloc(&db_graph, file.hdr.kmer_size,
                 file.hdr.num_of_cols, 0, kmers_in_hash);

  // In colour
  size_t words64_per_col = roundup_bits2words64(db_graph.ht.capacity);
  db_graph.node_in_cols = calloc2(words64_per_col * file.hdr.num_of_cols,
                                  sizeof(uint64_t));

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                             .boolean_covgs = false,
                             .must_exist_in_graph = false,
                             .must_exist_in_edges = NULL,
                             .empty_colours = false};

  graph_load(&file, gprefs, stats);

  if(add_pop_edges) status("Inferring edges from population...\n");
  else status("Inferring all missing edges...\n");

  // Read again
  fseek(file.fltr.fh, file.hdr_size, SEEK_SET);

  BinaryKmer bkmer;
  Edges edges[db_graph.num_of_cols];
  Covg covgs[db_graph.num_of_cols];

  size_t num_nodes_modified = 0;
  long edges_len = sizeof(Edges) * file.hdr.num_of_cols;

  while(graph_file_read(&file, &bkmer, covgs, edges))
  {
    if((add_all_edges && infer_all_edges(bkmer, edges, covgs, &db_graph)) ||
       (!add_all_edges && infer_pop_edges(bkmer, edges, covgs, &db_graph)))
    {
      if(fseek(file.fltr.fh, -edges_len, SEEK_CUR) != 0)
        die("fseek error: %s", path);
      fwrite(edges, 1, (size_t)edges_len, file.fltr.fh);
      num_nodes_modified++;
    }
  }

  char modified_str[100], kmers_str[100];
  ulong_to_str(num_nodes_modified, modified_str);
  ulong_to_str(db_graph.ht.unique_kmers, kmers_str);
  status("%s of %s (%.2f%%) nodes modified\n", modified_str, kmers_str,
         (100.0 * num_nodes_modified) / db_graph.ht.unique_kmers);

  seq_loading_stats_free(stats);
  free(db_graph.node_in_cols);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&file);

  return EXIT_SUCCESS;
}
