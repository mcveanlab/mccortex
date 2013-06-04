#include "global.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "add_read_paths.h"
#include "binary_format.h"

static const char usage[] =
"usage: ctx_add_paths <in.ctx> <mem> [OPTIONS]\n"
"  Thread reads through the graph.  Generates <in>.ctp file\n"
"  Options:\n"
"    --se_list <col> <in.list>\n"
"    --pe_list <col> <pe.list1> <pe.list2>\n";

int main(int argc, char* argv[])
{
  if(argc < 3) print_usage(usage, NULL);

  char *input_ctx_path = argv[1];
  char *mem_arg = argv[2];

  size_t mem_to_use = 0;

  // Check arguments
  if(!test_file_readable(input_ctx_path))
    print_usage(usage, "Cannot read input file: %s", input_ctx_path);

  if(!mem_to_integer(mem_arg, &mem_to_use) || mem_to_use == 0)
    print_usage(usage, "Invalid memory argument: %s", mem_arg);

  unsigned int col;
  int argi;
  for(argi = 3; argi < argc; argi++) {
    if(strcmp(argv[argi], "--se_list") == 0)
    {
      if(argi+2 >= argc)
        print_usage(usage, "--se_list <col> <input.falist> missing args");

      if(!parse_entire_uint(argv[argi+1], &col))
        print_usage(usage, "--se_list <col> <input.falist> invalid colour");

      check_colour_or_ctx_list(argv[argi+2], 0, false, true, 0);
      argi += 2;
    }
    else if(strcmp(argv[argi], "--pe_list") == 0)
    {
      if(argi+3 >= argc)
        print_usage(usage, "--pe_list <col> <in1.list> <in2.list> missing args");

      if(!parse_entire_uint(argv[argi+1], &col))
        print_usage(usage, "--pe_list <col> <in1.list> <in2.list> invalid colour");

      uint32_t num_files1, num_files2;
      num_files1 = check_colour_or_ctx_list(argv[argi+2], 0, false, true, 0);
      num_files2 = check_colour_or_ctx_list(argv[argi+3], 0, false, true, 0);
      if(num_files1 != num_files2)
        die("list mismatch [%s; %s]", argv[argi+2], argv[argi+3]);
      argi += 3;
    }
    else print_usage(usage, "Unknown argument: %s", argv[argi]);
  }

  // Probe binary to get kmer_size
  boolean is_binary = false;
  uint32_t kmer_size, num_of_cols;
  uint64_t num_kmers;

  if(!binary_probe(input_ctx_path, &is_binary, &kmer_size, &num_of_cols, &num_kmers))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  // Decide on memory
  size_t hash_kmers = num_kmers*(1/IDEAL_OCCUPANCY);

  size_t hash_memory = hash_kmers * sizeof(BinaryKmer) +
                       hash_kmers * sizeof(uint8_t*) +
                       hash_kmers * sizeof(Edges) +
                       hash_kmers * sizeof(uint8_t) +
                       hash_kmers * num_of_cols / 64;

  if(hash_memory > mem_to_use) {
    print_usage(usage, "Not enough memory; hash table requires %zu", hash_memory);
  }

  size_t path_memory = mem_to_use - hash_memory;

  message("Using %zu bytes hash; %zu bytes for paths\n", hash_memory, path_memory);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, hash_kmers);

  uint8_t *tmp = calloc(hash_kmers, sizeof(uint8_t)*2);
  if(tmp == NULL) die("Out of memory");

  db_graph.edges = tmp;
  db_graph.status = tmp + hash_kmers;

  size_t i, words64_per_col = round_bits_to_words64(hash_kmers);
  uint64_t *bkmer_cols = calloc(words64_per_col*NUM_OF_COLOURS, sizeof(uint64_t));
  if(bkmer_cols == NULL) die("Out of memory");

  uint64_t *ptr;
  for(ptr = bkmer_cols, i = 0; i < NUM_OF_COLOURS; i++, ptr += words64_per_col)
    db_graph.bkmer_in_cols[i] = ptr;

  db_graph.kmer_paths = malloc(hash_kmers * sizeof(uint64_t));
  if(db_graph.kmer_paths == NULL) die("Out of memory");
  memset(db_graph.kmer_paths, 0xff, path_memory);

  db_graph.kmer_paths = calloc(hash_kmers, sizeof(uint64_t));
  if(db_graph.kmer_paths == NULL) die("Out of memory");

  uint8_t *path_store = malloc(path_memory);
  if(path_store == NULL) die("Out of memory");

  binary_paths_init(&db_graph.pdata, path_store, path_memory);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  binary_load(input_ctx_path, &db_graph, 0, -1, true, false, stats);

  SeqLoadingPrefs prefs = {.into_colour = 0,
                           .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = false, .must_exist_in_colour = -1,
                           .empty_colours = false, .load_as_union = false,
                           .update_ginfo = true, .db_graph = &db_graph};

  // Parse input sequence
  for(argi = 3; argi < argc; argi++) {
    if(strcmp(argv[argi], "--se_list") == 0) {
      parse_entire_uint(argv[argi+1], &col);
      add_read_paths_to_graph(argv[argi+2], NULL, NULL, col, NULL, 0, prefs);
      argi += 2;
    }
    else if(strcmp(argv[argi], "--pe_list") == 0) {
      parse_entire_uint(argv[argi+1], &col);
      add_read_paths_to_graph(NULL, argv[argi+2], argv[argi+3], col,
                              NULL, 0, prefs);
      argi += 3;
    }
  }

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);
}
