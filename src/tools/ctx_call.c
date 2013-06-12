#include "global.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "add_read_paths.h"
#include "binary_format.h"
#include "graph_walker.h"
#include "shaded_caller.h"

static const char usage[] =
"usage: ctx_call <in.ctx> <mem> <out.bubbles.gz> [OPTIONS]\n"
"  Thread reads through the graph, calls bubbles.\n"
"  Options:\n"
"    --se_list <col> <in.list>\n"
"    --pe_list <col> <pe.list1> <pe.list2>\n";

#define NUM_THREADS 2

/*
void test_traverse(dBGraph *db_graph, char *str)
{
  BinaryKmer bkmer, bkey;
  binary_kmer_from_str(str, db_graph->kmer_size, bkmer);
  db_node_get_key(bkmer, db_graph->kmer_size, bkey);

  hkey_t node = hash_table_find(&db_graph->ht, bkey);
  Orientation or = db_node_get_orientation(bkmer, bkey);

  GraphWalker gwlk;
  graph_walker_alloc(&gwlk);
  graph_walker_init(&gwlk, db_graph, 0, node, or);

  printf("== test_traverse %s [dir:%s]\n", str, or == forward ? "fw" : "rv");

  char tmp[100];

  do {
    binary_kmer_to_str(gwlk.bkmer, db_graph->kmer_size, str);
    binary_kmer_to_str(db_node_bkmer(db_graph,gwlk.node), db_graph->kmer_size, tmp);
    printf("%s %s:%i\n", str, tmp, gwlk.orient);
  }
  while(graph_traverse(&gwlk));
  printf("\n");

  graph_walker_finish(&gwlk);

  graph_walker_dealloc(&gwlk);
}

void test(dBGraph *db_graph)
{
  Nucleotide bases[100] = {0};
  uint8_t bytes[100] = {0};

  bases[0] = 1;
  bases[1] = 1;
  bases[2] = 2;
  bases[3] = 3;
  bases[4] = 2;
  bases[5] = 3;

  // 3210 0032
  // 11100100 00001110
  // 228, 14

  size_t i;
  for(i = 0; i < 6; i++) printf("%c", binary_nuc_to_char(bases[i]));
  printf("\n");

  pack_bases(bytes, bases, 6);
  printf(" {%u,%u}\n", (uint32_t)bytes[0], (uint32_t)bytes[1]);
  unpack_bases(bytes, bases, 6);

  for(i = 0; i < 8; i++) printf("%c", binary_nuc_to_char(bases[i]));
  printf("\n");

  // TGCGTCGGCG
  // TCCGTCGGTG
  printf(" SEQ: TGCGTCGGCG\n");
  printf(" SEQ: TCCGTCGGTG\n");

  char tmp1[100], tmp2[100], tmp3[100], tmp4[100];
  strcpy(tmp1, "TGCGT");
  strcpy(tmp2, "TCCGT");
  strcpy(tmp3, "CGCCG");
  strcpy(tmp4, "CACCG");

  test_traverse(db_graph, tmp1);
  test_traverse(db_graph, tmp2);
  test_traverse(db_graph, tmp3);
  test_traverse(db_graph, tmp4);
}
*/

int main(int argc, char* argv[])
{
  if(argc < 7) print_usage(usage, NULL);

  char *input_ctx_path = argv[1];
  char *mem_arg = argv[2];
  char *out_path = argv[3];

  size_t mem_to_use = 0;

  // Check arguments
  if(!test_file_readable(input_ctx_path))
    print_usage(usage, "Cannot read input file: %s", input_ctx_path);

  if(!mem_to_integer(mem_arg, &mem_to_use) || mem_to_use == 0)
    print_usage(usage, "Invalid memory argument: %s", mem_arg);

  if(!test_file_writable(out_path))
    print_usage(usage, "Cannot write output file: %s", out_path);

  unsigned int col;
  int argi;
  for(argi = 4; argi < argc; argi++) {
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
  size_t req_num_kmers = num_kmers*(1.0/IDEAL_OCCUPANCY);
  size_t i, hash_kmers;
  size_t hash_mem = hash_table_mem(req_num_kmers, &hash_kmers);

  size_t graph_mem = hash_mem +
                     hash_kmers * sizeof(Edges) + // edges
                     hash_kmers * sizeof(uint64_t) * 2 + // kmer_paths
                     round_bits_to_bytes(hash_kmers) * num_of_cols + // in col
                     round_bits_to_bytes(hash_kmers) * 2; // visited fw/rv

  size_t thread_mem = round_bits_to_bytes(hash_kmers) * 2 * NUM_THREADS;

  if(graph_mem+thread_mem > mem_to_use) {
    print_usage(usage, "Not enough memory; hash table: %zu; threads: %zu",
                graph_mem, thread_mem);
  }

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, hash_kmers);

  size_t path_memory = mem_to_use - graph_mem - thread_mem;
  message("Using %zu bytes hash; %zu bytes for paths\n", graph_mem, path_memory);

  // Edges
  db_graph.edges = calloc(hash_kmers, sizeof(uint8_t));
  if(db_graph.edges == NULL) die("Out of memory");

  // In colour
  size_t words64_per_col = round_bits_to_words64(hash_kmers);
  uint64_t *bkmer_cols = calloc(words64_per_col*NUM_OF_COLOURS, sizeof(uint64_t));
  if(bkmer_cols == NULL) die("Out of memory");

  uint64_t *ptr;
  for(ptr = bkmer_cols, i = 0; i < NUM_OF_COLOURS; i++, ptr += words64_per_col)
    db_graph.node_in_cols[i] = ptr;

  // Visited
  // db_graph.visited = calloc(hash_kmers, sizeof(uint64_t));
  // if(db_graph.visited == NULL) die("Out of memory");

  // Paths
  db_graph.kmer_paths = malloc(hash_kmers * sizeof(uint64_t) * 2);
  if(db_graph.kmer_paths == NULL) die("Out of memory");
  memset(db_graph.kmer_paths, 0xff, hash_kmers * sizeof(uint64_t) * 2);

  uint8_t *path_store = malloc(path_memory);
  if(path_store == NULL) die("Out of memory");
  binary_paths_init(&db_graph.pdata, path_store, path_memory);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  binary_load(input_ctx_path, &db_graph, 0, -1, true, false, stats);

  //
  // Pack/Unpack test
  //
  // const char seq[] = "ATGGCGATAAGG";
  // char new_seq[100];

  // Nucleotide bases[100];
  // uint8_t compressed[100];

  // binary_nuc_from_str(bases, seq, strlen(seq));

  // pack_bases(compressed, bases, strlen(seq));
  // unpack_bases(compressed, bases, strlen(seq));

  // printf("From: %s\n", seq);
  // binary_nuc_to_str(bases, new_seq, strlen(seq));
  // printf("To  : %s\n", new_seq);

  // exit(-1);

  //
  //
  //

  SeqLoadingPrefs prefs = {.into_colour = 0,
                           .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = false, .must_exist_in_colour = -1,
                           .empty_colours = false, .load_as_union = false,
                           .update_ginfo = true, .db_graph = &db_graph};

  // Parse input sequence
  #define NUM_PASSES 1
  size_t rep;
  for(rep = 0; rep < NUM_PASSES; rep++)
  {
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
  }

  // db_graph_dump_paths_by_kmer(&db_graph);

  // Now call variants
  invoke_shaded_bubble_caller(&db_graph, out_path, NUM_THREADS);

  // test(&db_graph);

  free(db_graph.edges);
  free(bkmer_cols);
  // free(db_graph.visited);
  free(path_store);
  free(db_graph.kmer_paths);

  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  message("Done.\n");
}
