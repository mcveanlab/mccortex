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
"   -m <mem>     Memory to use\n"
"   -h <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   --merge      Merge corresponding colours from each binary\n"
"   --flatten    Dump into a single colour binary\n"
"\n"
"   --intersect <a.ctx>\n"
"     Only load the kmers that are in graph A.ctx. Can be specified multiple times.\n"
"     <a.ctx> is NOT merged into the output file.\n"
"\n"
"  Files can be specified with specific colours: samples.ctx:2,3\n"
"  Offset specifies where to load the first colour (merge only).\n";

static inline void remove_non_intersect_nodes(hkey_t node, Covg *covgs,
                                              Covg num, HashTable *ht)
{
  if(covgs[node] != num)
    hash_table_delete(ht, node);
}

int ctx_join(CmdArgs *args)
{
  cmd_accept_options(args, "mh");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  char *out_ctx_path;
  boolean merge = false, flatten = false;

  int argi;
  size_t num_intersect = 0;

  for(argi = 0; argi < argc; argi++) {
    if(strcasecmp(argv[argi],"--merge") == 0) {
      if(merge) warn("merge specified twice");
      merge = true;
    }
    else if(strcasecmp(argv[argi],"--flatten") == 0) {
      if(flatten) warn("flatten specified twice");
      flatten = true;
    }
    else if(strcasecmp(argv[argi],"--intersect") == 0)
    {
      if(argi+1 >= argc)
        print_usage(usage, "--intersect <A.ctx> requires and argument");
      num_intersect++;
      argi++;
    }
    else if(argv[argi][0] == '-') {
      print_usage(usage, "Unknown argument '%s'", argv[argi]);
    }
    else break;
  }

  // Store intersection files
  char *intersect_paths[num_intersect];
  int argj;
  num_intersect = 0;

  for(argj = 0; argj < argi; argj++) {
    if(strcasecmp(argv[argj],"--intersect") == 0) {
      intersect_paths[num_intersect++] = argv[argj+1];
      argj++;
    }
  }

  if(argc - argi < 2)
    print_usage(usage, "Please specify output and input binaries");

  out_ctx_path = argv[argi++];

  // argi .. argend-1 are binaries to load
  uint32_t num_binaries = argc - argi;
  char **binary_paths = argv + argi;

  // Check all binaries are valid binaries with matching kmer size
  boolean is_binary = false;
  uint32_t i, kmer_size, kmer_size2;
  uint64_t ctx_num_kmers;
  uint32_t ctx_max_cols[num_binaries], ctx_num_cols[num_binaries];
  char *path;

  for(i = 0; i < num_binaries; i++)
  {
    // Strip off offset (e.g. 12:in.ctx)
    for(path = binary_paths[i]; *path >= '0' && *path <= '9'; path++) {}

    if(path > binary_paths[i] && *path == ':') path++;
    else path = binary_paths[i];

    if(!binary_probe(path, &is_binary, &kmer_size2, &ctx_num_cols[i],
                     &ctx_max_cols[i], &ctx_num_kmers)) {
      print_usage(usage, "Cannot read input binary file: %s", binary_paths[i]);
    }
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", binary_paths[i]);

    if(i == 0)
      kmer_size = kmer_size2;
    else if(kmer_size != kmer_size2)
      print_usage(usage, "Kmer sizes don't match [%u vs %u]", kmer_size, kmer_size2);
  }

  // Probe intersection files
  uint32_t intersect_num_cols, intersect_max_cols;
  uint64_t intersect_num_kmers, min_intersect_num_kmers;

  for(i = 0; i < num_intersect; i++)
  {
    if(!binary_probe(intersect_paths[i], &is_binary, &kmer_size2, &intersect_num_cols,
                     &intersect_max_cols, &intersect_num_kmers)) {
      print_usage(usage, "Cannot read intersect binary file: %s", intersect_paths[i]);
    }
    else if(!is_binary)
      print_usage(usage, "Intersect binary file isn't valid: %s", intersect_paths[i]);
  
    if(i == 0) min_intersect_num_kmers = intersect_num_kmers;
    else if(intersect_num_kmers < min_intersect_num_kmers)
    {
      // Put smallest intersection binary first
      char *tmpstr;
      SWAP(intersect_paths[i], intersect_paths[0], tmpstr);
      min_intersect_num_kmers = intersect_num_kmers;
    }
  }

  // Pick hash table size
  size_t kmers_in_hash;
  size_t mem_per_kmer = sizeof(BinaryKmer) + sizeof(Covg) + sizeof(Edges);

  if(num_intersect > 0)
  {
    size_t ideal_occupancy = min_intersect_num_kmers / IDEAL_OCCUPANCY;
    if(args->num_kmers_set) {
      char tmp[50];
      if(args->num_kmers < min_intersect_num_kmers) {
        die("Need at least -h %s", bytes_to_str(min_intersect_num_kmers,1,tmp));
      } else if(args->num_kmers < ideal_occupancy) {
        warn("Less than ideal hash occupancy (-h %zu < %zu [70%%])",
             args->num_kmers, ideal_occupancy);
      }
      kmers_in_hash = args->num_kmers;
    }
    else kmers_in_hash = ideal_occupancy;

    hash_table_mem(kmers_in_hash, &kmers_in_hash);

    size_t mem_needed = kmers_in_hash * mem_per_kmer;
    char mem_needed_str[50], mem_set_str[50];
    bytes_to_str(mem_needed, 1, mem_needed_str);
    bytes_to_str(args->mem_to_use, 1, mem_set_str);

    if(args->mem_to_use_set && mem_needed > args->mem_to_use)
      die("Requires %s memory, you set -m %s", mem_needed_str, mem_set_str);

    message("[memory] graph: %s\n", mem_needed_str);
  }
  else
    kmers_in_hash = cmd_get_kmers_in_hash(args, mem_per_kmer);

  message("\n");

  // Check out_ctx_path is writable
  if(!test_file_writable(out_ctx_path))
    print_usage(usage, "Cannot write to output: %s", out_ctx_path);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, kmers_in_hash);
  db_graph.col_edges = calloc(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc(db_graph.ht.capacity, sizeof(Covg));

  // Load intersection binaries
  if(num_intersect > 0)
  {
    SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                             .boolean_covgs = true,
                             .load_seq = true,
                             .quality_cutoff = 0, .ascii_fq_offset = 0,
                             .homopolymer_cutoff = 0,
                             .remove_dups_se = false, .remove_dups_pe = false,
                             .load_binaries = true,
                             .must_exist_in_graph = false, .empty_colours = false,
                             .update_ginfo = false,
                             .db_graph = &db_graph};

    binary_load(intersect_paths[i], &db_graph, &prefs, NULL, NULL);

    if(num_intersect > 1)
    {
      prefs.must_exist_in_graph = true;

      for(i = 1; i < num_intersect; i++)
        binary_load(intersect_paths[i], &db_graph, &prefs, NULL, NULL);

      // Remove nodes where covg != num_intersect
      HASH_TRAVERSE(&db_graph.ht, remove_non_intersect_nodes, db_graph.col_covgs,
                    num_intersect, &db_graph.ht);
    }

    // Zero covgs and edges
    memset(db_graph.col_edges, 0, db_graph.ht.capacity * sizeof(Edges));
    memset(db_graph.col_covgs, 0, db_graph.ht.capacity * sizeof(Covg));

    message("Loaded intersection set\n\n");
  }

  binaries_merge(out_ctx_path, binary_paths, num_binaries,
                 ctx_num_cols, ctx_max_cols,
                 merge, flatten, num_intersect > 0, &db_graph);

  free(db_graph.col_edges);
  free(db_graph.col_covgs);
  db_graph_dealloc(&db_graph);
  message("Done.\n");
  return EXIT_SUCCESS;
}

//
/*

  message("Using kmer size %u; Creating %u colour binary\n\n",
          kmer_size, output_colours);

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
       .boolean_covgs = false,
       .load_seq = true,
       .quality_cutoff = 0, .ascii_fq_offset = 0,
       .homopolymer_cutoff = 0,
       .remove_dups_se = false, .remove_dups_pe = false,
       .load_binaries = true,
       .must_exist_in_graph = false,
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
    uint32_t load_colours[num_binaries][max_cols];
    for(i = 0; i < num_binaries; i++)
      binary_parse_colour_array(binary_paths[i], load_colours[i], ctx_max_cols[i]);

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
      {
        graph_info_merge(output_header.ginfo + output_colour,
                         tmpheader.ginfo + load_colours[i][j]);
      }
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
  */
