
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
"usage: "CMD" join [-m <mem>] <out.ctx> <in1.ctx> [in2.ctx ...]\n"
"  Merge cortex binaries.\n"
"\n"
" Options:\n"
"   --merge    Merge corresponding colours from each binary\n"
"   --flatten  Dump into a single colour binary\n"
"\n"
" Files can be specified with specific colours: samples.ctx:2,3\n";

static inline void dump_empty_bkmer(hkey_t node, dBGraph *db_graph,
                                    char *buf, size_t mem, FILE *fh)
{
  fwrite(db_node_bkmer(db_graph, node), sizeof(BinaryKmer), 1, fh);
  fwrite(buf, 1, mem, fh);
}

static void dump_empty_binary(dBGraph *db_graph, FILE *fh, uint32_t num_of_cols)
{
  size_t mem = num_of_cols * (sizeof(Covg)+sizeof(Edges));
  char buf[mem];
  memset(buf, 0, mem);
  HASH_TRAVERSE(&db_graph->ht, dump_empty_bkmer, db_graph, buf, mem, fh);
}

int ctx_join(CmdArgs *args)
{
  cmd_accept_options(args, "m");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  uint64_t mem_to_use = args->mem_to_use;

  char *out_ctx_path;
  boolean merge = false, flatten = false;

  int i, argstart;

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
  int num_binaries = argc - argstart;

  // Check all binaries are valid binaries with matching kmer size
  boolean is_binary = false;
  uint32_t kmer_size, kmer_size2, bin_num_cols[num_binaries];
  uint32_t max_cols = 0, sum_cols = 0;
  uint64_t num_kmers;
  char *ctx_path;

  for(i = 0; i < num_binaries; i++)
  {
    ctx_path = argv[argstart+i];

    if(!binary_probe(ctx_path, &is_binary, &kmer_size2, bin_num_cols+i, &num_kmers))
      print_usage(usage, "Cannot read input binary file: %s", ctx_path);
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", ctx_path);

    if(i == 0)
      kmer_size = kmer_size2;
    else if(kmer_size != kmer_size2)
      print_usage(usage, "Kmer sizes don't match [%u vs %u]", kmer_size, kmer_size2);

    max_cols = MAX2(bin_num_cols[i], max_cols);
    sum_cols += bin_num_cols[i];
  }

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  uint32_t output_colours;

  if(flatten) output_colours = 1;
  else if(merge) output_colours = max_cols;
  else output_colours = sum_cols;

  // Pick hash table size
  size_t mem_per_kmer, kmers_in_hash, hash_mem, graph_mem;

  mem_per_kmer = sizeof(BinaryKmer) + sizeof(Covg) + sizeof(Edges);
  hash_mem = hash_table_mem2(mem_to_use / mem_per_kmer, &kmers_in_hash);
  graph_mem = kmers_in_hash * mem_per_kmer;

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, max_cols, kmers_in_hash);

  db_graph.col_edges = calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc(db_graph.ht.capacity, sizeof(Covg));

  // Print mem usage
  char graph_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  message("[memory]  graph: %s\n", graph_mem_str);
  hash_table_print_stats(&db_graph.ht);

  message("Using kmer size %u; Creating %u colour binary\n",
          kmer_size, output_colours);

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                           .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_colour = -1,
                           .empty_colours = num_binaries == 1,
                           .update_ginfo = false,
                           .db_graph = &db_graph};
  //

  if(output_colours == 1)
  {
    // e.g. flatten
    for(i = 0; i < num_binaries; i++)
    {
      prefs.into_colour = 0;
      binary_load(argv[argstart+i], &db_graph, &prefs, stats);
    }
    
    binary_dump_graph(out_ctx_path, &db_graph, CURR_CTX_VERSION, NULL, 0, 1);
  }
  else
  {
    // Construct binary header
    BinaryFileHeader header = {.version = CURR_CTX_VERSION,
                               .kmer_size = db_graph.kmer_size,
                               .num_of_bitfields = NUM_BITFIELDS_IN_BKMER,
                               .num_of_cols = output_colours,
                               .num_of_kmers = db_graph.ht.unique_kmers};

    header.ginfo = malloc(sizeof(GraphInfo) * output_colours);

    Colour j, output_colour = 0;
    for(j = 0; j < output_colours; j++) graph_info_alloc(header.ginfo + j);

    for(i = 0; i < num_binaries; i++)
    {
      prefs.into_colour = 0;
      binary_load(argv[argstart+i], &db_graph, &prefs, stats);

      // if merge set output_colour to zero
      output_colour *= !merge;

      for(j = 0; j < bin_num_cols[i]; j++, output_colour++) {
        graph_info_merge(header.ginfo + output_colour, db_graph.ginfo + j);
      }
    }

    db_graph_set_cols(&db_graph, 1);

    FILE *fh = fopen(out_ctx_path, "w");
    if(fh == NULL) die("Cannot open output ctx file: %s", out_ctx_path);

    size_t header_size = binary_write_header(fh, &header);

    // print file outline
    dump_empty_binary(&db_graph, fh, output_colours);

    if(merge)
    {
      for(output_colour = 0; output_colour < output_colours; output_colour++)
      {
        prefs.into_colour = 0;
        memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
        memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));
      
        for(i = 0; i < num_binaries; i++)
        {
          if(output_colour < bin_num_cols[i])
          {
            binary_load_colour(argv[argstart+i], &db_graph, &prefs, stats,
                               output_colour);
          }
        }
        fseek(fh, header_size, SEEK_SET);
        binary_dump_colour(&db_graph, 0, output_colour, output_colours, fh);
      }
    }
    else
    {
      output_colour = 0;
      for(i = 0; i < num_binaries; i++)
      {
        prefs.into_colour = 0;
        memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
        memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));

        for(j = 0; j < bin_num_cols[i]; j++)
        {
          binary_load_colour(argv[argstart+i], &db_graph, &prefs, stats, j);

          fseek(fh, header_size, SEEK_SET);
          binary_dump_colour(&db_graph, 0, output_colour, output_colours, fh);
          output_colour++;
        }
      }
    }

    fclose(fh);
    free(header.ginfo);
  }

  hash_table_print_stats(&db_graph.ht);

  message("Dumped %zu kmers in %u colour%s\n", (size_t)db_graph.ht.unique_kmers,
          output_colours, output_colours != 1 ? "s" : "");
  message("Done.\n");

  seq_loading_stats_free(stats);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
