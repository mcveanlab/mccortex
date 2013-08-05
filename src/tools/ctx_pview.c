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

static const char usage[] = "usage: "CMD" pview [options] <in.ctp>\n";

int ctx_pview(CmdArgs *args)
{
  cmd_accept_options(args, "m");
  // cmd_require_options(args, "m", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 1) print_usage(usage, NULL);

  char *input_paths_file = argv[0];

  // probe paths file to get kmer size
  boolean valid_paths_file = false;
  PathFileHeader pheader = {.capacity = 0};

  if(!paths_file_probe(input_paths_file, &valid_paths_file, &pheader))
  {
    print_usage(usage, "Cannot find .ctp file: %s", input_paths_file);
  }

  if(!valid_paths_file)
    die("Invalid .ctp file: %s", input_paths_file);

  char num_paths_str[100], path_bytes_str[100], kmers_with_paths_str[100];
  ulong_to_str(pheader.num_of_paths, num_paths_str);
  bytes_to_str(pheader.num_path_bytes, 1, path_bytes_str);
  ulong_to_str(pheader.num_kmers_with_paths, kmers_with_paths_str);

  // Print header
  message("version: %u\n", pheader.version);
  message("kmer size: %u\n", pheader.kmer_size);
  message("colours: %u\n", pheader.num_of_cols);
  message("paths: %s\n", num_paths_str);
  message("bytes: %s\n", path_bytes_str);
  message("kmers: %s\n", kmers_with_paths_str);

  uint32_t col;
  for(col = 0; col < pheader.num_of_cols; col++)
    message(" colour %u: %s\n", col, pheader.sample_names[col].buff);

  // Decide on memory
  size_t kmers_in_hash, ideal_capacity, req_num_kmers;
  size_t hash_mem, graph_mem, path_mem;

  ideal_capacity = pheader.num_kmers_with_paths / IDEAL_OCCUPANCY;
  req_num_kmers = args->num_kmers_set ? args->num_kmers : ideal_capacity;
  hash_mem = hash_table_mem(req_num_kmers, &kmers_in_hash);

  graph_mem = hash_mem +
                     kmers_in_hash * sizeof(uint64_t); // kmer_paths

  path_mem = args->mem_to_use - graph_mem;

  char graph_mem_str[100], path_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(path_mem, 1, path_mem_str);

  message("[memory]  graph: %s;  paths: %s\n", graph_mem_str, path_mem_str);

  if(kmers_in_hash < pheader.num_kmers_with_paths) {
    print_usage(usage, "Not enough kmers in the hash, require: %s "
                       "(set bigger -h <kmers> or -m <mem>)", kmers_with_paths_str);
  }
  else if(kmers_in_hash < ideal_capacity)
    warn("Low memory for binary size (require: %s)", kmers_with_paths_str);

  if(args->mem_to_use_set && graph_mem > args->mem_to_use)
    die("Not enough memory (please increase -m <mem>)");

  // Allocate memory
  // db graph is required to store the end position for each kmer list
  dBGraph db_graph;
  db_graph_alloc(&db_graph, pheader.kmer_size, 1, kmers_in_hash);

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  binary_paths_init(&db_graph.pdata, path_store, path_mem, pheader.num_of_cols);

  // Pretend we've read all the kmers in
  db_graph.num_of_cols_used = pheader.num_of_cols;

  // Add kmers as reading
  boolean add_kmers = true;

  paths_format_read(input_paths_file, &pheader, &db_graph,
                    &db_graph.pdata, add_kmers);

  db_graph_dump_paths_by_kmer(&db_graph);

  free(path_store);
  free((void *)db_graph.kmer_paths);

  paths_header_dealloc(&pheader);
  db_graph_dealloc(&db_graph);

  message("Done.\n");
  return EXIT_SUCCESS;
}
