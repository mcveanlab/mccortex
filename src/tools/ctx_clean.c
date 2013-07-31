#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "binary_format.h"
#include "supernode.h"

// DEV: add option to dump csv of supernode covg before cleaning

// clean cortex binaries
static const char usage[] =
"usage: "CMD" clean [options] <out.ctx> <in.ctx> [in2.ctx ...]\n"
"  Clean a cortex binary.\n"
"  Options:\n"
"    -m <mem>        Memory to use\n"
"    -h <hash-size>  Kmers in the hash table (e.g. 1G ~ 1 billion)\n"
"    -g <genome>     Genome size\n"
"    --threshold <T> Cleaning threshold\n"
"    --tips          Clip tips\n"
"    --supernodes    Remove low coverage supernode\n"
"                    (requires -g <G> or --threshold <T>)\n"
"  With no options only supernode cleaning is done\n";

static void printsupernodes(hkey_t node, dBGraph *db_graph, uint64_t *visited)
{
  uint32_t kmer_size = db_graph->kmer_size;
  if(!bitset_has(visited, node))
  {
    bitset_set(visited, node);
    char bkmer[MAX_KMER_SIZE+1];
    binary_kmer_to_str(db_node_bkmer(db_graph, node), kmer_size, bkmer);
    printf("%s\n", bkmer);

    Supernode snode;
    supernode_load(node, db_graph, &snode);
    char bkmer1[MAX_KMER_SIZE+1], bkmer2[MAX_KMER_SIZE+1];
    ConstBinaryKmerPtr bkmerptr1 = db_node_bkmer(db_graph, snode.start_node);
    ConstBinaryKmerPtr bkmerptr2 = db_node_bkmer(db_graph, snode.end_node);
    binary_kmer_to_str(bkmerptr1, kmer_size, bkmer1);
    binary_kmer_to_str(bkmerptr2, kmer_size, bkmer2);
    printf("%s:%i -> %s:%i\n", bkmer1, snode.start_orient,
                               bkmer2, snode.end_orient);
  }
}

int ctx_clean(CmdArgs *args)
{
  cmd_accept_options(args, "mhg");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  // Check cmdline args
  boolean tip_cleaning = false, supernode_cleaning = false;
  uint32_t threshold = 0;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++) {
    if(strcmp(argv[argi],"--tips") == 0) tip_cleaning = true;
    else if(strcmp(argv[argi],"--supernodes") == 0) supernode_cleaning = true;
    else if(strcmp(argv[argi],"--threshold") == 0) {
      if(argi + 1 >= argc || !parse_entire_uint(argv[argi+1], &threshold) ||
         threshold <= 1) {
        print_usage(usage, "--threshold <T> needs an integer argument > 1");
      }
      argi++;
    }
    else print_usage(usage, "Unknown argument: %s", argv[argi]);
  }

  if(argc - argi < 2) print_usage(usage, "Please give file names");

  char *out_ctx_path = argv[argi++];
  char **binary_paths = argv + argi;
  uint32_t i, j, num_binaries = argc - argi;

  // default settings
  if(!tip_cleaning && !supernode_cleaning)
    supernode_cleaning = true;

  if(supernode_cleaning && threshold == 0 && !args->genome_size_set)
    print_usage(usage, "supernode cleaning requires --threshold <T> or "
                       "-g <genome-size>");

  if(args->genome_size_set && threshold > 0)
    print_usage(usage, "Only one one --threshold <T> and -g <genome-size> "
                       "needed: <T> is calculated from <genome-size>");

  // Probe binary files
  boolean is_binary = false;
  uint32_t kmer_size, kmer_size2;
  uint32_t num_of_cols = 0, max_cols = 0, ctx_max_col, ctx_num_cols[num_binaries];
  uint64_t num_kmers;

  for(i = 0; i < num_binaries; i++)
  {
    if(!binary_probe(binary_paths[i], &is_binary, &kmer_size2, &ctx_num_cols[i],
                     &ctx_max_col, &num_kmers)) {
      print_usage(usage, "Cannot read input binary file: %s", binary_paths[i]);
    }
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", binary_paths[i]);

    if(i == 0)
      kmer_size = kmer_size2;
    else if(kmer_size != kmer_size2)
      print_usage(usage, "Kmer sizes don't match [%u vs %u]", kmer_size, kmer_size2);

    num_of_cols += ctx_num_cols[i];
    max_cols = MAX2(max_cols, ctx_num_cols[i]);

    printf("%s has %u colours\n", binary_paths[i], ctx_num_cols[i]);
  }

  uint32_t load_colours[num_binaries][max_cols];
  for(i = 0; i < num_binaries; i++)
    binary_parse_colour_array(binary_paths[i], load_colours[i], ctx_max_col);

  // Pick hash table size
  size_t mem_per_kmer = sizeof(BinaryKmer) + sizeof(Covg) + sizeof(Edges);
  size_t kmers_in_hash = cmd_get_kmers_in_hash(args, mem_per_kmer);

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, kmers_in_hash);
  db_graph.col_edges = calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc(db_graph.ht.capacity, sizeof(Covg));

  // Load binary into a single colour
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs
    = {.into_colour = 0, .merge_colours = true,
       .boolean_covgs = false,
       .load_seq = true,
       .quality_cutoff = 0, .ascii_fq_offset = 0,
       .homopolymer_cutoff = 0,
       .remove_dups_se = false, .remove_dups_pe = false,
       .load_binaries = true,
       .must_exist_in_graph = false,
       .empty_colours = false,
       .update_ginfo = true,
       .db_graph = &db_graph};

  // Construct cleaned binary header
  BinaryFileHeader tmpheader;
  BinaryFileHeader output_header = {.version = CURR_CTX_VERSION,
                                    .kmer_size = db_graph.kmer_size,
                                    .num_of_bitfields = NUM_BITFIELDS_IN_BKMER,
                                    .num_of_cols = num_of_cols,
                                    .num_of_kmers = db_graph.ht.unique_kmers};

  binary_header_alloc(&tmpheader, max_cols);
  binary_header_alloc(&output_header, num_of_cols);

  uint32_t output_colour = 0;
  for(i = 0; i < num_binaries; i++) {
    binary_load(binary_paths[i], &db_graph, &prefs, stats, &tmpheader);
    for(j = 0; j < ctx_num_cols[i]; j++, output_colour++)
      graph_info_merge(output_header.ginfo + output_colour, tmpheader.ginfo + j);
  }

  hash_table_print_stats(&db_graph.ht);
  uint64_t *visited = malloc(round_bits_to_words64(db_graph.ht.capacity) * sizeof(uint64_t));

  // DEV: Clean
  HASH_TRAVERSE(&db_graph.ht, printsupernodes, &db_graph, visited);

  FILE *fh = fopen(out_ctx_path, "w");
    if(fh == NULL) die("Cannot open output ctx file: %s", out_ctx_path);

  size_t header_size = binary_write_header(fh, &output_header);

  // Free header resources
  binary_header_dealloc(&tmpheader);
  binary_header_dealloc(&output_header);

  dump_empty_binary(&db_graph, fh, num_of_cols);

  // load, clean and dump graph one colour at a time
  prefs.must_exist_in_graph = true;

  output_colour = 0;
  for(i = 0; i < num_binaries; i++)
  {
    for(j = 0; j < num_of_cols; j++)
    {
      memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
      memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));
      binary_load_colour(binary_paths[i], &db_graph, &prefs, stats,
                         load_colours[i][j]);
      fseek(fh, header_size, SEEK_SET);
      binary_dump_colour(&db_graph, 0, output_colour, num_of_cols, fh);
      output_colour++;
    }
  }

  fclose(fh);
 
  message("Dumped %zu kmers in %u colour%s into: %s\n",
          (size_t)db_graph.ht.unique_kmers, num_of_cols,
          num_of_cols != 1 ? "s" : "", out_ctx_path);
  message("Done.\n");

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
