#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_walker.h"
#include "supernode.h"
#include "repeat_walker.h"

static const char usage[] =
"usage: "CMD" contigs [options] <input.ctx>\n"
"  Pull out contigs, print statistics\n"
"  Options: [ -m <mem> | -n <kmers> | -p <paths.ctp> ]\n"
"           [ --ncontigs <N> | --print | --colour <c> ]\n";

int ctx_contigs(CmdArgs *args)
{
  cmd_accept_options(args, "mnpo", usage);
  int argc = args->argc;
  char **argv = args->argv;

  // char str[100];
  // printf("Test NaN: %s INF: %s\n", num_to_str(NAN, 2, str),
  //                                  num_to_str(INFINITY, 2, str));

  size_t ncontigs = 10000, colour = 0;
  boolean print_contigs = false;

  while(argc > 0 && argv[0][0] == '-') {
    if(strcmp(argv[0],"--ncontigs") == 0) {
      unsigned long tmp;
      if(argc == 1 || !parse_entire_ulong(argv[1], &tmp))
        print_usage(usage, "--ncontigs <N> requires an integer argument");
      ncontigs = tmp;
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--colour") == 0 || strcmp(argv[0],"--color") == 0) {
      unsigned long tmp;
      if(argc == 1 || !parse_entire_ulong(argv[1], &tmp))
        print_usage(usage, "--colour <c> requires an integer argument");
      colour = tmp;
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--print") == 0) {
      print_contigs = true;
      argv++; argc--;
    }
    else print_usage(usage, "Unknown argument: %s", argv[0]);
  }

  if(argc != 1) print_usage(usage, NULL);
  char *input_ctx_path = argv[0];

  // Seed random
  seed_random();

  // probe binary
  GraphFileReader file = INIT_GRAPH_READER;
  int ret = graph_file_open(&file, input_ctx_path, false);

  if(ret == 0)
    print_usage(usage, "Cannot read input graph file: %s", input_ctx_path);
  else if(ret < 0)
    print_usage(usage, "Input graph file isn't valid: %s", input_ctx_path);

  // probe paths files
  boolean valid_paths_file = false;
  PathFileHeader pheader = INIT_PATH_FILE_HDR;
  size_t i;

  if(args->num_ctp_files == 0)
    status("No path files (.ctp) to load");
  else if(args->num_ctp_files > 1)
    print_usage(usage, "Cannot load >1 .ctp file at the moment [use pmerge]");

  for(i = 0; i < args->num_ctp_files; i++)
  {
    if(!paths_file_probe(args->ctp_files[i], &valid_paths_file, &pheader))
      print_usage(usage, "Cannot read .ctp file: %s", args->ctp_files[i]);
    else if(!valid_paths_file)
      die("Invalid .ctp file: %s", args->ctp_files[i]);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem;

  bits_per_kmer = sizeof(Edges)*8 + file.hdr.num_of_cols + sizeof(uint64_t)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        file.hdr.num_of_kmers, false);

  graph_mem = hash_table_mem(kmers_in_hash, false, NULL) +
              (bits_per_kmer * kmers_in_hash) / 8;
  path_mem = pheader.num_path_bytes;

  char path_mem_str[100];
  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  if(graph_mem + path_mem > args->mem_to_use) die("Not enough memory");

  //
  // Output file if printing
  //
  FILE *fout = stdout;
  if(args->output_file_set) {
    if(!print_contigs)
      warn("Ignoring --out <out> argument (maybe you forgot --print ?)");
    else if((fout = fopen(args->output_file, "r")) == NULL)
      die("Cannot open output file: %s", args->output_file);
  }

  dBGraph db_graph;
  GraphWalker wlk;

  db_graph_alloc(&db_graph, file.hdr.kmer_size, file.hdr.num_of_cols, 1, kmers_in_hash);
  graph_walker_alloc(&wlk);

  size_t node_bit_fields = round_bits_to_words64(db_graph.ht.capacity);

  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.node_in_cols = calloc2(node_bit_fields * file.hdr.num_of_cols,
                                  sizeof(uint64_t));
  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(uint64_t));

  RepeatWalker rptwlk;
  walker_alloc(&rptwlk, db_graph.ht.capacity, 22); // 4MB

  uint8_t *path_store = malloc2(path_mem);
  path_store_init(&db_graph.pdata, path_store, path_mem, file.hdr.num_of_cols);

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = true};

  graph_load(&file, &prefs, stats);
  seq_loading_stats_free(stats);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  for(i = 0; i < args->num_ctp_files; i++) {
    paths_format_read(args->ctp_files[i], &pheader, &db_graph,
                      &db_graph.pdata, false);
  }

  status("Traversing graph...\n");

  hkey_t node;
  Orientation orient;
  Nucleotide lost_nuc;

  size_t njunc;
  size_t total_len = 0, total_junc = 0, contigs_outdegree[5] = {0}, nloop = 0;
  size_t *lengths = malloc2(ncontigs * sizeof(size_t));
  size_t *junctions = malloc2(ncontigs * sizeof(size_t));
  double density, max_density = 0;
  size_t max_len = 0, max_junc = 0;

  size_t len, path_cap = 4096;
  dBNode *nodes = malloc2(path_cap * sizeof(dBNode));

  for(i = 0; i < ncontigs; i++)
  {
    do { node = db_graph_rand_node(&db_graph); }
    while(!db_node_has_col(&db_graph, node, colour));

    nodes[0].key = node;
    nodes[0].orient = FORWARD;
    len = 1;
    njunc = 0;

    for(orient = 0; orient < 2; orient++)
    {
      if(orient == 1) {
        supernode_reverse(nodes, len);
        node = nodes[len-1].key;
      }

      graph_walker_init(&wlk, &db_graph, colour, colour, node, orient);
      lost_nuc = binary_kmer_first_nuc(wlk.bkmer, db_graph.kmer_size);

      while(graph_traverse(&wlk))
      {
        if(walker_attempt_traverse(&rptwlk, &wlk, wlk.node, wlk.orient, wlk.bkmer))
        {
          graph_walker_node_add_counter_paths(&wlk, lost_nuc);
          lost_nuc = binary_kmer_first_nuc(wlk.bkmer, db_graph.kmer_size);
          if(len == path_cap) {
            path_cap *= 2;
            nodes = realloc2(nodes, path_cap * sizeof(dBNode));
          }
          nodes[len].key = wlk.node;
          nodes[len].orient = wlk.orient;
          len++;
        }
        else break;
      }

      njunc += wlk.fork_count;
      nloop += rptwlk.nbloom_entries;
      graph_walker_finish(&wlk);
      walker_fast_clear(&rptwlk, nodes, len);
    }

    if(print_contigs) {
      fprintf(fout, ">contig%zu\n", i);
      db_nodes_print(nodes, len, &db_graph, fout);
      putc('\n', fout);
    }

    size_t outdegree = edges_get_outdegree(db_graph.col_edges[wlk.node], wlk.orient);
    contigs_outdegree[outdegree]++;
    lengths[i] = len;
    junctions[i] = njunc;
    total_len += len;
    total_junc += njunc;

    density = (double)njunc / len;
    max_density = MAX2(max_density, density);
    max_len = MAX2(max_len, len);
    max_junc = MAX2(max_junc, njunc);
  }

  if(args->output_file_set && print_contigs)
    fclose(fout);

  status("\n");

  qsort(lengths, ncontigs, sizeof(size_t), cmp_size);
  qsort(junctions, ncontigs, sizeof(size_t), cmp_size);

  double median_len = MEDIAN(lengths, ncontigs);
  double median_junc = MEDIAN(junctions, ncontigs);

  char ncontigs_str[90];
  char len_mean_str[90], len_median_str[90], len_max_str[90], len_total_str[90];
  char jnc_mean_str[90], jnc_median_str[90], jnc_max_str[90], jnc_total_str[90];

  num_to_str((double)total_len / ncontigs, 1, len_mean_str);
  num_to_str((double)total_junc / ncontigs, 1, jnc_mean_str);
  num_to_str(median_len, 1, len_median_str);
  num_to_str(median_junc, 1, jnc_median_str);
  num_to_str(max_len, 1, len_max_str);
  num_to_str(max_junc, 1, jnc_max_str);
  num_to_str(total_len, 1, len_total_str);
  num_to_str(total_junc, 1, jnc_total_str);

  num_to_str(ncontigs, 1, ncontigs_str);

  status("Pulled out %s contigs", ncontigs_str);
  status("Lengths: mean: %s, median: %s, max: %s, total: %s",
         len_mean_str, len_median_str, len_max_str, len_total_str);
  status("Junctions: mean: %s, median: %s, max: %s, total: %s",
         jnc_mean_str, jnc_median_str, jnc_max_str, jnc_total_str);
  status("Max junction density: %.2f\n", max_density);
  status("Contigs looping back to a kmer: %zu [%.2f%%]\n", nloop,
         (100.0 * nloop) / ncontigs);

  timestamp(ctx_msg_out);
  message(" Contig outdegree: ");

  for(i = 0; i <= 4; i++) {
    message(" %zu: %zu [%.2f%%]", i, contigs_outdegree[i],
            (100.0*contigs_outdegree[i])/(2*ncontigs));
  }
  message("\n");

  free(lengths);
  free(junctions);
  free(nodes);

  // free(visited);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free((void*)db_graph.kmer_paths);
  free(path_store);

  walker_dealloc(&rptwlk);
  paths_header_dealloc(&pheader);
  graph_walker_dealloc(&wlk);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&file);

  return EXIT_SUCCESS;
}
