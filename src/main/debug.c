#include "global.h"
#include "file_util.h"
#include "db_graph.h"
#include "gpath_hash.h"
#include "gpath_reader.h"
#include "gpath_save.h"
#include "gpath_checks.h"

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

  // Create a path store that tracks path counts
  gpath_store_alloc(&db_graph.gpstore,
                    db_graph.num_of_cols, db_graph.ht.capacity,
                    ONE_MEGABYTE, true);

  // Create path hash table for fast lookup
  GPathHash phash;
  gpath_hash_alloc(&phash, &db_graph.gpstore, ONE_MEGABYTE);

  // DEV: check sample names match

  // Load path files
  gpath_reader_load(&pfile, false, &db_graph);
  gpath_reader_close(&pfile);

  hash_table_print_stats(&db_graph.ht);

  // Write output file
  gpath_save(gzout, out_path, &db_graph);
  gzclose(gzout);

  // Checks
  // gpath_checks_all_paths(&db_graph);
  gpath_checks_counts(&db_graph);

  // Clean up
  db_graph_dealloc(&db_graph);
  cortex_destroy();

  return EXIT_SUCCESS;
}
