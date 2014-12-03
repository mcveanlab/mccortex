#include "global.h"
#include "cmd.h"
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
  cmd_init(argc, argv);

  if(argc != 3) die("usage: ./debug <in.ctp> <in.ctx>");

  const char *out_path = argv[2];

  GPathReader pfile;
  memset(&pfile, 0, sizeof(GPathReader));
  gpath_reader_open(&pfile, argv[1], true);
  status("Got file with %zu colours", pfile.ncolours);

  size_t i, kmer_size = 7, ncols = 3;

  gpath_reader_check(&pfile, kmer_size, ncols);
  gzFile gzout = futil_gzopen_create(out_path, "w");

  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, ncols, 1, 1024, DBG_ALLOC_EDGES);

  // Create a path store that tracks path counts
  gpath_store_alloc(&db_graph.gpstore,
                    db_graph.num_of_cols, db_graph.ht.capacity,
                    ONE_MEGABYTE, true, false);

  // Create path hash table for fast lookup
  gpath_hash_alloc(&db_graph.gphash, &db_graph.gpstore, ONE_MEGABYTE);

  // Set sample names
  for(i = 0; i < pfile.ncolours; i++) {
    const char *sample_name = gpath_reader_get_sample_name(&pfile, i);
    ctx_assert(sample_name != NULL);
    strbuf_set(&db_graph.ginfo[i].sample_name, sample_name);
  }

  // Load path files, add kmers that are missing
  gpath_reader_load(&pfile, GPATH_ADD_MISSING_KMERS, &db_graph);

  hash_table_print_stats(&db_graph.ht);

  // Write output file
  gpath_save(gzout, out_path, 1, true, &pfile.json, 1, &db_graph);
  gzclose(gzout);

  // Checks
  // gpath_checks_all_paths(&db_graph, 2); // use two threads
  gpath_checks_counts(&db_graph);

  // Clean up
  gpath_reader_close(&pfile);
  db_graph_dealloc(&db_graph);
  cortex_destroy();

  return EXIT_SUCCESS;
}
