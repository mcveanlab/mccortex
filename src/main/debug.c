#include "global.h"
#include "file_util.h"
#include "db_graph.h"
#include "gpath_reader.h"
#include "gpath_save.h"

int main(int argc, char **argv)
{
  (void)argc; (void)argv;
  cortex_init();

  if(argc != 3) die("usage: ./debug <in.ctp> <out.ctp>");

  const char *out_path = argv[2];

  GPathReader pfile;
  memset(&pfile, 0, sizeof(GPathReader));
  gpath_reader_open(&pfile, argv[1], true);

  size_t kmer_size = 7, ncols = 3;

  gpath_reader_check(&pfile, kmer_size, ncols);
  gzFile gzout = futil_gzopen_output(out_path);

  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, ncols, 1, 1024);
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));

  gpath_store_alloc(&db_graph.gpstore, 10000,
                    db_graph.ht.capacity, db_graph.num_of_cols,
                    true, true);

  // DEV: check sample names match

  gpath_reader_load(&pfile, false, &db_graph);
  gpath_reader_close(&pfile);

  hash_table_print_stats(&db_graph.ht);

/*  char *ctx_path = argv[1];
  GraphFileReader file = INIT_GRAPH_READER;
  graph_file_open(&file, ctx_path, true); // true => errors are fatal

  GraphLoadingPrefs prefs = {.db_graph = &db_graph,
                             .boolean_covgs = false,
                             .must_exist_in_graph = false,
                             .must_exist_in_edges = NULL,
                             .empty_colours = true};

  graph_load(&file, prefs, NULL);

  hkey_t nodes[1];
  BinaryKmer bkmer = binary_kmer_from_str("ATATATATCTAGATATATATCTATATATAAA", kmer_size);
  BinaryKmer bkey = binary_kmer_get_key(bkmer, kmer_size);
  nodes[0] = hash_table_find(&db_graph.ht, bkey);
  ctx_assert(nodes[0] != HASH_NOT_FOUND);

*/

  gpath_save(gzout, out_path, &db_graph);
  gzclose(gzout);

  db_graph_dealloc(&db_graph);

  cortex_destroy();

  return EXIT_SUCCESS;
}
