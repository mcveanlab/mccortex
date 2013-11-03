#include "global.h"
#include <time.h>
#include <pthread.h>

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_store.h"
#include "path_format.h"

static const char usage[] =
"usage: "CMD" pmerge [-m <mem> | -h <kmer> | -f <in.ctx>] <out.ctp> <in1.ctp> ...\n"
"  Merge paths.  One of either -h <kmers> or -f <in.ctx> is required to specify\n"
"  the size of the hash table\n";

int ctx_pmerge(CmdArgs *args)
{
  cmd_accept_options(args, "mhf");
  cmd_accept_options(args, "m");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  if(!args->num_kmers_set && !args->file_set)
    print_usage(usage, "Please specify -h <num-kmers> or -f <in.ctx>");
  if(args->num_kmers_set && args->file_set)
    print_usage(usage, "Please specify only ONE of -h <num-kmers> or -f <in.ctx>");

  char *input_ctx_path = argv[0];
  char *out_ctp_path = argv[1];

  // probe binary
  uint64_t num_kmers = args->num_kmers;
  uint32_t ctp_num_of_cols, ctp_kmer_size;
  GraphFileReader file = INIT_GRAPH_READER;

  if(args->file_set)
  {
    int ret = graph_file_open(&file, input_ctx_path, false);

    if(ret == 0)
      print_usage(usage, "Cannot read input graph file: %s", input_ctx_path);
    else if(ret < 0)
      print_usage(usage, "Input graph file isn't valid: %s", input_ctx_path);

    num_kmers = file.hdr.num_of_kmers;
  }

  uint64_t max_ctp_path_bytes = 0;
  PathFileHeader pheader = INIT_PATH_FILE_HDR;

  int argi;
  for(argi = 2; argi < argc; argi++)
  {
    boolean valid_paths_file = false;
    const char *input_paths_file = argv[argi];

    if(!file_exists(input_paths_file))
    {
      print_usage(usage, "Cannot find ctp file: %s", input_paths_file);
    }
    else if(!paths_file_probe(input_paths_file, &valid_paths_file, &pheader))
    {
      print_usage(usage, "Cannot read .ctp file: %s", input_paths_file);
    }
    else if(!valid_paths_file)
      die("Invalid .ctp file: %s", input_paths_file);
    else if(num_kmers < pheader.num_kmers_with_paths) {
      if(args->file_set) die("ctp file has more kmers than hash table!");
      else {
        print_usage(usage, "Please set a larger -h <kmers> (needs to be > %zu)",
                    (size_t)pheader.num_kmers_with_paths);
      }
    }

    // If we didn't load a binary, assume num of colours from first .ctp file
    if(argi == 2) {
      ctp_kmer_size = pheader.kmer_size;
      ctp_num_of_cols = pheader.num_of_cols;
    }
    else
    {
      if(pheader.num_of_cols != ctp_num_of_cols) {
        die("Number of colours in .ctp files does not match: %s [%u] vs %s [%u]",
            argv[2], ctp_num_of_cols, argv[argi], pheader.num_of_cols);
      } else if(pheader.kmer_size != ctp_kmer_size) {
        die("Kmer size in .ctp files does not match: %s [%u] vs %s [%u]",
            argv[2], ctp_kmer_size, argv[argi], pheader.kmer_size);
      }
    }
    max_ctp_path_bytes = MAX2(max_ctp_path_bytes, pheader.num_path_bytes);
  }

  if(args->file_set) {
    if(ctp_num_of_cols != file.hdr.num_of_cols) {
      warn("Number of colours in .ctp files does not match .ctx [%u vs %u]",
           ctp_num_of_cols, file.hdr.num_of_cols);
    } else if(ctp_kmer_size != file.hdr.kmer_size) {
      warn("Kmer size in .ctp files does not match .ctx [%u vs %u]",
           ctp_kmer_size, file.hdr.kmer_size);
    }
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, path_mem;
  char path_mem_str[100];

  bits_per_kmer = sizeof(uint64_t)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        file.hdr.num_of_kmers, true);
  path_mem = args->mem_to_use -
             kmers_in_hash * (sizeof(BinaryKmer)+sizeof(uint64_t));

  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  // set up tmp space
  PathStore tmp_pdata;
  uint8_t *tmp_path_store = malloc2(max_ctp_path_bytes);
  path_store_init(&tmp_pdata, tmp_path_store, max_ctp_path_bytes, file.hdr.num_of_cols);

  // Set up graph and PathStore
  dBGraph db_graph;

  db_graph_alloc(&db_graph, ctp_kmer_size, ctp_num_of_cols, 0, kmers_in_hash);

  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  path_store_init(&db_graph.pdata, path_store, path_mem, file.hdr.num_of_cols);

  //
  // Set up file header
  //
  pheader.version = CTX_PATH_FILEFORMAT;
  pheader.kmer_size = db_graph.kmer_size;
  pheader.num_of_cols = file.hdr.num_of_cols;
  pheader.capacity = 0;

  paths_header_alloc(&pheader, file.hdr.num_of_cols);

  // Recursively load/add paths
  paths_format_read(argv[2], &pheader, &db_graph, &db_graph.pdata, true);

  for(argi = 3; argi < argc; argi++) {
    paths_format_merge(argv[argi], &pheader, &db_graph,
                       &db_graph.pdata, &tmp_pdata, true);
  }

  // Dump paths file
  FILE *fout = fopen(out_ctp_path, "w");
  if(fout == NULL) die("Cannot open output file: %s", out_ctp_path);
  setvbuf(fout, NULL, _IOFBF, CTX_BUF_SIZE);
  paths_header_update(&pheader, &db_graph.pdata);
  paths_format_write_header(&pheader, fout);
  paths_format_write_optimised_paths(&db_graph, fout);
  fclose(fout);

  free((void *)db_graph.kmer_paths);
  free(path_store);

  paths_header_dealloc(&pheader);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&file);

  status("Paths written to: %s\n", out_ctp_path);

  return EXIT_SUCCESS;
}
