#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "add_read_paths.h"
#include "graph_format.h"
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
  PathFileHeader pheader = INIT_PATH_FILE_HDR;

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
  printf("version: %u\n", pheader.version);
  printf("kmer size: %u\n", pheader.kmer_size);
  printf("colours: %u\n", pheader.num_of_cols);
  printf("paths: %s\n", num_paths_str);
  printf("bytes: %s\n", path_bytes_str);
  printf("kmers starting paths: %s\n", kmers_with_paths_str);

  uint32_t col;
  for(col = 0; col < pheader.num_of_cols; col++) {
    printf(" colour %u: %s\n", col, pheader.sample_names[col].buff);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;

  bits_per_kmer = sizeof(uint64_t) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        pheader.num_kmers_with_paths, true);

  // Allocate memory
  // db graph is required to store the end position for each kmer list
  dBGraph db_graph;
  db_graph_alloc(&db_graph, pheader.kmer_size, 1, 0, kmers_in_hash);

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  uint8_t *path_store = malloc2(pheader.num_path_bytes);
  path_store_init(&db_graph.pdata, path_store,
                  pheader.num_path_bytes, pheader.num_of_cols);

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

  return EXIT_SUCCESS;
}
