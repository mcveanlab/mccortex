#include "global.h"
#include <time.h>
#include <pthread.h>

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"
#include "binary_paths.h"
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
  uint32_t ctp_num_of_cols, ctx_num_of_cols, ctx_max_col;
  uint32_t ctp_kmer_size, ctx_kmer_size;

  if(args->file_set)
  {
    boolean is_binary = false;

    if(!binary_probe(input_ctx_path, &is_binary, &ctx_kmer_size,
                     &ctx_num_of_cols, &ctx_max_col, &num_kmers)) {
      print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
    } else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);
  }

  uint64_t max_ctp_path_bytes = 0;

  int argi;
  for(argi = 2; argi < argc; argi++)
  {
    boolean valid_paths_file = false;
    uint64_t ctp_num_paths, ctp_num_path_bytes, ctp_num_path_kmers;
    uint32_t tmp_kmer_size, tmp_num_of_cols;
    const char *input_paths_file = argv[argi];

    if(!file_exists(input_paths_file))
    {
      print_usage(usage, "Cannot find ctp file: %s", input_paths_file);
    }
    else if(!paths_format_probe(input_paths_file, &valid_paths_file,
                                &tmp_kmer_size, &tmp_num_of_cols, &ctp_num_paths,
                                &ctp_num_path_bytes, &ctp_num_path_kmers))
    {
      print_usage(usage, "Cannot read .ctp file: %s", input_paths_file);
    }
    else if(!valid_paths_file)
      die("Invalid .ctp file: %s", input_paths_file);
    else if(num_kmers < ctp_num_path_kmers) {
      if(args->file_set) die("ctp file has more kmers than hash table!");
      else {
        print_usage(usage, "Please set a larger -h <kmers> (needs to be > %zu)",
                    (size_t)ctp_num_path_kmers);
      }
    }

    // If we didn't load a binary, assume num of colours from first .ctp file
    if(argi == 2) {
      ctp_kmer_size = tmp_kmer_size;
      ctp_num_of_cols = tmp_num_of_cols;
    }
    else
    {
      if(tmp_num_of_cols != ctp_num_of_cols) {
        die("Number of colours in .ctp files does not match: %s [%u] vs %s [%u]",
            argv[2], ctp_num_of_cols, argv[argi], tmp_num_of_cols);
      } else if(tmp_kmer_size != ctp_kmer_size) {
        die("Kmer size in .ctp files does not match: %s [%u] vs %s [%u]",
            argv[2], ctp_kmer_size, argv[argi], tmp_kmer_size);
      }
    }
    max_ctp_path_bytes = MAX2(max_ctp_path_bytes, ctp_num_path_bytes);
  }

  if(args->file_set) {
    if(ctp_num_of_cols != ctx_num_of_cols) {
      warn("Number of colours in .ctp files does not match .ctx [%u vs %u]",
           ctp_num_of_cols, ctx_num_of_cols);
    } else if(ctp_kmer_size != ctx_kmer_size) {
      warn("Kmer size in .ctp files does not match .ctx [%u vs %u]",
           ctp_kmer_size, ctx_kmer_size);
    }
  }

  // Decide on memory
  size_t kmers_in_hash, req_num_kmers = num_kmers*(1.0/IDEAL_OCCUPANCY);
  size_t hash_mem = hash_table_mem(req_num_kmers, &kmers_in_hash);
  size_t path_mem = args->mem_to_use - hash_mem;

  if(args->mem_to_use_set && hash_mem+path_mem > args->mem_to_use)
    die("Not enough memory (please increase -m <mem>)");

  char hash_mem_str[100], path_mem_str[100], num_kmers_str[100];
  ulong_to_str(req_num_kmers, num_kmers_str);
  bytes_to_str(hash_mem, 1, hash_mem_str);
  bytes_to_str(path_mem, 1, path_mem_str);

  message("[memory]  graph: %s;  paths: %s\n", hash_mem_str, path_mem_str);

  // set up tmp space
  PathStore tmp_pdata;
  uint8_t *tmp_path_store = malloc(max_ctp_path_bytes);
  if(tmp_path_store == NULL) die("Out of memory");
  binary_paths_init(&tmp_pdata, tmp_path_store, max_ctp_path_bytes, ctx_num_of_cols);

  // Set up graph and PathStore
  dBGraph db_graph;

  db_graph_alloc(&db_graph, ctp_kmer_size, ctp_num_of_cols, req_num_kmers);

  db_graph.kmer_paths = malloc(db_graph.ht.capacity * sizeof(uint64_t));
  if(db_graph.kmer_paths == NULL) die("Out of memory");
  memset((void*)db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(uint64_t));

  uint8_t *path_store = malloc(path_mem);
  if(path_store == NULL) die("Out of memory");
  binary_paths_init(&db_graph.pdata, path_store, path_mem, ctx_num_of_cols);

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
