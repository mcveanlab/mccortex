#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "add_read_paths.h"
#include "binary_format.h"
#include "path_format.h"
#include "graph_walker.h"
#include "bubble_caller.h"

static const char usage[] = "usage: "CMD" pview [options] <in.ctx>\n";

int ctx_pview(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 1) print_usage(usage, NULL);

  size_t mem_to_use = args->mem_to_use;
  if(!args->mem_to_use_set) print_usage(usage, "-m <M> required");

  char *input_ctx_path = argv[0];

  // probe paths file to get kmer size
  char *input_paths_file = malloc(strlen(input_ctx_path)+4);
  paths_format_filename(input_ctx_path, input_paths_file);
  boolean valid_paths_file = false;
  uint64_t ctp_num_paths = 0, ctp_num_path_bytes = 0, ctp_num_path_kmers = 0;
  uint32_t ctp_kmer_size = 0, ctp_num_of_cols = 0;

  if(!paths_format_probe(input_paths_file, &valid_paths_file,
                         &ctp_kmer_size, &ctp_num_of_cols, &ctp_num_paths,
                         &ctp_num_path_bytes, &ctp_num_path_kmers))
  {
    print_usage(usage, "Cannot find .ctp file: %s", input_paths_file);
  }

  if(!valid_paths_file)
    die("Invalid .ctp file: %s", input_paths_file);

  // Decide on memory
  size_t req_num_kmers = ctp_num_path_kmers*(1.0/IDEAL_OCCUPANCY);
  size_t hash_kmers;
  size_t hash_mem = hash_table_mem(req_num_kmers, &hash_kmers);

  size_t graph_mem = hash_mem +
                     hash_kmers * sizeof(uint64_t); // kmer_paths

  // Allocate memory
  // db graph is required to store the end position for each kmer list
  dBGraph db_graph;
  db_graph_alloc(&db_graph, ctp_kmer_size, 1, hash_kmers);

  size_t path_mem = mem_to_use - graph_mem;

  char graph_mem_str[100], path_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(path_mem, 1, path_mem_str);

  message("[memory]  graph: %s;  paths: %s\n", graph_mem_str, path_mem_str);

  // Paths
  db_graph.kmer_paths = malloc(hash_kmers * sizeof(uint64_t));
  if(db_graph.kmer_paths == NULL) die("Out of memory");
  memset((void*)db_graph.kmer_paths, 0xff, hash_kmers * sizeof(uint64_t));

  uint8_t *path_store = malloc(path_mem);
  if(path_store == NULL) die("Out of memory");
  binary_paths_init(&db_graph.pdata, path_store, path_mem, ctp_num_of_cols);

  // Pretend we've read all the kmers in
  db_graph.num_of_cols_used = ctp_num_of_cols;

  // Add kmers as reading
  paths_format_read(&db_graph, &db_graph.pdata, NULL, true, input_paths_file);

  db_graph_dump_paths_by_kmer(&db_graph);

  free(path_store);
  free((void *)db_graph.kmer_paths);

  db_graph_dealloc(&db_graph);
  free(input_paths_file);

  message("Done.\n");
  return EXIT_SUCCESS;
}
