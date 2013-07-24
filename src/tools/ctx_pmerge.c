#include "global.h"
#include <time.h>
#include <pthread.h>

#include "cmd.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"
#include "binary_paths.h"
#include "path_format.h"

static const char usage[] =
"usage: "CMD" pmerge [options] <in.ctx> <out.ctp> <in1.ctp> ...\n";

int ctx_pmerge(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  char *input_ctx_path = argv[0];
  char *out_ctp_path = argv[1];

  // probe binary
  boolean is_binary = false;
  uint32_t kmer_size, num_of_cols;
  uint64_t num_kmers;

  if(!binary_probe(input_ctx_path, &is_binary, &kmer_size, &num_of_cols, &num_kmers))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  uint64_t max_ctp_path_bytes = 0;

  int argi;
  for(argi = 2; argi < argc; argi++)
  {
    boolean valid_paths_file = false;
    uint64_t ctp_num_paths, ctp_num_path_bytes, ctp_num_path_kmers;
    uint32_t ctp_kmer_size, ctp_num_of_cols;
    const char *input_paths_file = argv[argi];

    if(!file_exists(input_paths_file))
    {
      print_usage(usage, "Cannot find ctp file: %s", input_paths_file);
    }
    else if(!paths_format_probe(input_paths_file, &valid_paths_file,
                                &ctp_kmer_size, &ctp_num_of_cols, &ctp_num_paths,
                                &ctp_num_path_bytes, &ctp_num_path_kmers))
    {
      print_usage(usage, "Cannot read .ctp file: %s", input_paths_file);
    }
    else if(!valid_paths_file)
      die("Invalid .ctp file: %s", input_paths_file);
    else if(ctp_num_of_cols != num_of_cols)
      die("Number of colours in .ctp does not match .ctx");
    else if(ctp_kmer_size != kmer_size)
      die("Kmer size in .ctp does not match .ctx");
  
    max_ctp_path_bytes = MAX2(max_ctp_path_bytes, ctp_num_path_bytes);
  }

  // Decide on memory
  size_t hash_kmers, req_num_kmers = num_kmers*(1.0/IDEAL_OCCUPANCY);
  size_t hash_mem = hash_table_mem(req_num_kmers, &hash_kmers);
  size_t path_mem = args->mem_to_use - hash_mem;

  // set up tmp space
  PathStore tmp_pdata;
  uint8_t *tmp_path_store = malloc(max_ctp_path_bytes);
  if(tmp_path_store == NULL) die("Out of memory");
  binary_paths_init(&tmp_pdata, tmp_path_store, max_ctp_path_bytes, num_of_cols);

  // Set up graph and PathStore
  dBGraph db_graph;

  db_graph_alloc(&db_graph, kmer_size, num_of_cols, req_num_kmers);

  db_graph.kmer_paths = malloc(db_graph.ht.capacity * sizeof(uint64_t));
  if(db_graph.kmer_paths == NULL) die("Out of memory");
  memset((void*)db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(uint64_t));

  uint8_t *path_store = malloc(path_mem);
  if(path_store == NULL) die("Out of memory");
  binary_paths_init(&db_graph.pdata, path_store, path_mem, num_of_cols);

  // Recursively load/add paths
  paths_format_read(&db_graph, &db_graph.pdata, NULL, true, argv[2]);

  for(argi = 3; argi < argc; argi++)
    paths_format_read(&db_graph, &db_graph.pdata, &tmp_pdata, true, argv[argi]);

  // Dump paths file
  paths_format_write(&db_graph, &db_graph.pdata, out_ctp_path);

  free((void *)db_graph.kmer_paths);
  free(path_store);

  db_graph_dealloc(&db_graph);

  message("  Paths written to: %s\n", out_ctp_path);
  message("Done.");

  return EXIT_SUCCESS;
}
