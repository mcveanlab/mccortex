
#include "global.h"

#include "string_buffer.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "graph_info.h"
#include "db_node.h"
#include "binary_format.h"


// Given (A,B,C) are ctx binaries, A:1 means colour 1 in A,
// {A:1,B:0} is loading A:1 and B:0 into a single colour
//
// Default behaviour is to load colours consecutively
//   output: {A:0},{A:1},{B:0},{B:1},{B:2},{C:0},{C:1}
//
// --flatten
//   All colours into one
//   output: {A:0,A:1,B:0,B:1,B:2,C:0,C:1}
//
// --merge
//   Join input colours
//   output: {A:0,B:0,C:0},{A:1,B:1,C:1},{B:2}

static const char usage[] =
"usage: "CMD" join [options] <out.ctx> [offset:]in1.ctx[:1,2,4-5] [in2.ctx ...]\n"
"  Merge cortex binaries.  \n"
"\n"
"  Options:\n"
"   -m <mem>   Memory to use\n"
"   -h <kmers> Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   --merge    Merge corresponding colours from each binary\n"
"   --flatten  Dump into a single colour binary\n"
"\n"
"  Files can be specified with specific colours: samples.ctx:2,3\n"
"  Offset specifies where to load the first colour (merge only).\n";

int ctx_join(CmdArgs *args)
{
  cmd_accept_options(args, "mh");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  char *out_ctx_path;
  boolean merge = false, flatten = false;

  int argstart;

  for(argstart = 0; argstart < argc; argstart++) {
    if(strcasecmp(argv[argstart],"--merge") == 0) {
      if(merge) warn("merge specified twice");
      merge = true;
    }
    else if(strcasecmp(argv[argstart],"--flatten") == 0) {
      if(flatten) warn("flatten specified twice");
      flatten = true;
    }
    else if(argv[argstart][0] == '-') {
      print_usage(usage, "Unknown argument '%s'", argv[argstart]);
    }
    else break;
  }

  if(argc - argstart < 2)
    print_usage(usage, "Please specify output and input binaries");

  out_ctx_path = argv[argstart++];

  // argstart .. argend-1 are binaries to load
  uint32_t num_binaries = argc - argstart;
  char **binary_paths = argv + argstart;
  uint32_t i, offsets[num_binaries];
  uint32_t ctx_max_cols[num_binaries], ctx_num_cols[num_binaries];

  // Check all binaries are valid binaries with matching kmer size
  boolean is_binary = false;
  uint32_t kmer_size, kmer_size2;
  uint32_t max_cols = 0, sum_cols = 0;
  uint64_t num_kmers;
  char *ptr, *endptr;

  for(i = 0; i < num_binaries; i++)
  {
    for(ptr = binary_paths[i]; *ptr >= '0' && *ptr <= '9'; ptr++) {}

    if(ptr > binary_paths[i] && *ptr == ':') {
      if(!merge) {
        die("Cannot specify binary offsets without --merge option (%s)",
            binary_paths[i]);
      }
      offsets[i] = strtoul(binary_paths[i], &endptr, 10);
      binary_paths[i] = ptr+1;
    } else {
      offsets[i] = 0;
    }

    if(!binary_probe(binary_paths[i], &is_binary, &kmer_size2, &ctx_num_cols[i],
                     &ctx_max_cols[i], &num_kmers)) {
      print_usage(usage, "Cannot read input binary file: %s", binary_paths[i]);
    }
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", binary_paths[i]);

    if(i == 0)
      kmer_size = kmer_size2;
    else if(kmer_size != kmer_size2)
      print_usage(usage, "Kmer sizes don't match [%u vs %u]", kmer_size, kmer_size2);

    max_cols = MAX2(offsets[i] + ctx_num_cols[i], max_cols);
    sum_cols += ctx_num_cols[i];

    printf("%s has %u colours\n", binary_paths[i], ctx_num_cols[i]);
  }

  uint32_t output_colours;

  if(flatten) output_colours = 1;
  else if(merge) output_colours = max_cols;
  else output_colours = sum_cols;

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

  message("Using kmer size %u; Creating %u colour binary\n\n",
          kmer_size, output_colours);

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs
    = {.into_colour = 0, .merge_colours = true,
       .load_seq = true,
       .quality_cutoff = 0, .ascii_fq_offset = 0,
       .homopolymer_cutoff = 0,
       .remove_dups_se = false, .remove_dups_pe = false,
       .load_binaries = true,
       .must_exist_in_colour = -1,
       .empty_colours = num_binaries == 1 && output_colours == 1,
       .update_ginfo = false,
       .db_graph = &db_graph};
  //

  if(output_colours == 1)
  {
    // e.g. flatten

    for(i = 0; i < num_binaries; i++)
      binary_load(binary_paths[i], &db_graph, &prefs, stats, NULL);

    hash_table_print_stats(&db_graph.ht);    
    binary_dump_graph(out_ctx_path, &db_graph, CURR_CTX_VERSION, NULL, 0, 1);
  }
  else
  {
    // Construct binary header
    BinaryFileHeader tmpheader;
    BinaryFileHeader output_header = {.version = CURR_CTX_VERSION,
                                      .kmer_size = db_graph.kmer_size,
                                      .num_of_bitfields = NUM_BITFIELDS_IN_BKMER,
                                      .num_of_cols = output_colours,
                                      .num_of_kmers = db_graph.ht.unique_kmers};

    binary_header_alloc(&tmpheader, max_cols);
    binary_header_alloc(&output_header, output_colours);

    Colour j, output_colour = 0;
    for(i = 0; i < num_binaries; i++)
    {
      binary_load(binary_paths[i], &db_graph, &prefs, stats, &tmpheader);
      if(merge) output_colour = 0;
      for(j = 0; j < ctx_num_cols[i]; j++, output_colour++)
        graph_info_merge(output_header.ginfo + output_colour, tmpheader.ginfo + j);
    }

    FILE *fh = fopen(out_ctx_path, "w");
    if(fh == NULL) die("Cannot open output ctx file: %s", out_ctx_path);

    size_t header_size = binary_write_header(fh, &output_header);

    // Free header resources
    binary_header_dealloc(&tmpheader);
    binary_header_dealloc(&output_header);

    // print file outline
    message("Generated merged hash table\n\n");
    hash_table_print_stats(&db_graph.ht);
    dump_empty_binary(&db_graph, fh, output_colours);

    uint32_t load_colours[num_binaries][max_cols];
    for(i = 0; i < num_binaries; i++)
      binary_parse_colour_array(binary_paths[i], load_colours[i], ctx_max_cols[i]);

    if(merge)
    {
      for(output_colour = 0; output_colour < output_colours; output_colour++)
      {
        memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
        memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));
      
        boolean data_loaded_in_col = false;
        for(i = 0; i < num_binaries; i++)
        {
          uint32_t ctx_col;
          if(output_colour >= offsets[i] &&
             (ctx_col = output_colour - offsets[i]) < ctx_num_cols[i])
          {
            binary_load_colour(binary_paths[i], &db_graph, &prefs, stats,
                               load_colours[i][ctx_col]);
            data_loaded_in_col = true;
          }
        }
        if(data_loaded_in_col) {
          message("Dumping into colour %u...\n\n", output_colour);
          fseek(fh, header_size, SEEK_SET);
          binary_dump_colour(&db_graph, 0, output_colour, output_colours, fh);
        }
      }
    }
    else
    {
      output_colour = 0;
      for(i = 0; i < num_binaries; i++)
      {
        for(j = 0; j < ctx_num_cols[i]; j++, output_colour++)
        {
          memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
          memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));

          binary_load_colour(binary_paths[i], &db_graph, &prefs, stats,
                             load_colours[i][j]);

          message("Dumping into colour %u...\n\n", output_colour);
          fseek(fh, header_size, SEEK_SET);
          binary_dump_colour(&db_graph, 0, output_colour, output_colours, fh);
        }
      }
    }

    fclose(fh);
  }

  message("Dumped %zu kmers in %u colour%s into: %s\n",
          (size_t)db_graph.ht.unique_kmers, output_colours,
          output_colours != 1 ? "s" : "", out_ctx_path);
  message("Done.\n");

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
