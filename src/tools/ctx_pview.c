#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_walker.h"
#include "bubble_caller.h"

static const char usage[] = "usage: "CMD" pview [options] <in.ctp>\n";

int ctx_pview(CmdArgs *args)
{
  cmd_accept_options(args, "mn", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 1) print_usage(usage, NULL);

  char *input_paths_file = argv[0];

  // Open paths file
  PathFileReader pfile = INIT_PATH_READER;
  path_file_open(&pfile, input_paths_file, true);

  PathFileHeader *phdr = &pfile.hdr;
  char num_paths_str[100], path_bytes_str[100], kmers_with_paths_str[100];
  ulong_to_str(phdr->num_of_paths, num_paths_str);
  bytes_to_str(phdr->num_path_bytes, 1, path_bytes_str);
  ulong_to_str(phdr->num_kmers_with_paths, kmers_with_paths_str);

  // Print header
  printf("version: %u\n", phdr->version);
  printf("kmer size: %u\n", phdr->kmer_size);
  printf("colours: %u\n", phdr->num_of_cols);
  printf("paths: %s\n", num_paths_str);
  printf("bytes: %s\n", path_bytes_str);
  printf("kmers starting paths: %s\n", kmers_with_paths_str);

  size_t col;
  for(col = 0; col < phdr->num_of_cols; col++) {
    printf(" colour %zu: %s\n", col, phdr->sample_names[col].buff);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(uint64_t) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        phdr->num_kmers_with_paths,
                                        true, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  // Allocate memory
  // db graph is required to store the end position for each kmer list
  dBGraph db_graph;
  db_graph_alloc(&db_graph, phdr->kmer_size, phdr->num_of_cols, 0, kmers_in_hash);

  path_file_set_graph_sample_names(&pfile, &db_graph);

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(uint64_t));

  path_store_alloc(&db_graph.pdata, phdr->num_path_bytes, 0, phdr->num_of_cols);

  // Pretend we've read all the kmers in
  db_graph.num_of_cols_used = phdr->num_of_cols;

  // Add kmers as reading
  boolean add_kmers = true;

  paths_format_load(&pfile, &db_graph, add_kmers);
  db_graph_dump_paths_by_kmer(&db_graph);

  free((void *)db_graph.kmer_paths);

  path_store_dealloc(&db_graph.pdata);
  db_graph_dealloc(&db_graph);
  path_file_dealloc(&pfile);

  return EXIT_SUCCESS;
}
