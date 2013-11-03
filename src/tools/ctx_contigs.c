#include "global.h"
#include <sys/time.h> // srand

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

static const char usage[] =
"usage: "CMD" contigs [options] <input.ctx>\n"
"  Pull out contigs, print statistics\n"
"  Options: [ -m <mem> | -h <kmers> | -p <paths.ctp> ]\n"
"           [ --nsamples <N> | --print | --colour <c> ]\n";

int ctx_contigs(CmdArgs *args)
{
  cmd_accept_options(args, "mhp");
  int argc = args->argc;
  char **argv = args->argv;

  size_t num_samples = 10000, colour = 0;
  boolean print_contigs = false;

  while(argc > 0 && argv[0][0] == '-') {
    if(strcmp(argv[0],"--nsamples") == 0) {
      unsigned long tmp;
      if(argc == 1 || !parse_entire_ulong(argv[1], &tmp))
        print_usage(usage, "--nsamples <N> requires an integer argument");
      num_samples = tmp;
      argv += 2; argc -= 2;
    }
    else if(strcmp(argv[0],"--colour") == 0) {
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
  struct timeval time;
  gettimeofday(&time, NULL);
  srand((((time.tv_sec ^ getpid()) * 1000000) + time.tv_usec));

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

  // Get starting bkmer
  // BinaryKmer bkmer, bkey;
  // if(strlen(start_kmer) != kmer_size) die("length of kmer does not match kmer_size");
  // bkmer = binary_kmer_from_str(start_kmer, kmer_size);

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

  uint64_t *visited = calloc2(2 * node_bit_fields, sizeof(uint64_t));

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

  // Find start node
  // hkey_t node;
  // Orientation orient;
  // db_node_get_key(bkmer, kmer_size, bkey);
  // node = hash_table_find(&db_graph.ht, bkey);
  // orient = db_node_get_orientation(bkmer, bkey);
  Nucleotide lost_nuc;

  // char bkmerstr[MAX_KMER_SIZE+1];

  hkey_t node;
  Orientation orient;

  size_t j, njunc;
  size_t total_len = 0, total_junc = 0, dead_ends = 0, nloop = 0;
  size_t lengths[num_samples], junctions[num_samples];
  double density, max_density = 0;
  size_t max_len = 0, max_junc = 0;

  size_t len, path_cap = 4096;
  hkey_t *nodes = malloc2(path_cap * sizeof(hkey_t));
  Orientation *orients = malloc2(path_cap * sizeof(Orientation));

  for(i = 0; i < num_samples; i++)
  {
    node = db_graph_rand_node(&db_graph);
    nodes[0] = node;
    orients[0] = FORWARD;
    len = 1;
    njunc = 0;

    for(orient = 0; orient < 2; orient++)
    {
      if(orient == 1) {
        supernode_reverse(nodes, orients, len);
        node = nodes[len-1];
      }

      graph_walker_init(&wlk, &db_graph, colour, colour, node, orient);
      lost_nuc = binary_kmer_first_nuc(wlk.bkmer, db_graph.kmer_size);

      while(graph_traverse(&wlk))
      {
        if(db_node_has_traversed(visited, wlk.node, wlk.orient)){ nloop++; break;}
        db_node_set_traversed(visited, wlk.node, wlk.orient);
        graph_walker_node_add_counter_paths(&wlk, lost_nuc);
        lost_nuc = binary_kmer_first_nuc(wlk.bkmer, db_graph.kmer_size);
        if(len == path_cap) {
          path_cap *= 2;
          nodes = realloc2(nodes, path_cap * sizeof(hkey_t));
          orients = realloc2(orients, path_cap * sizeof(Orientation));
        }
        nodes[len] = wlk.node;
        orients[len] = wlk.orient;
        len++;
      }

      njunc += wlk.fork_count;
      graph_walker_finish(&wlk);

      for(j = 0; j < len; j++) db_node_fast_clear_traversed(visited, nodes[j]);
    }

    if(print_contigs) {
      fprintf(stdout, ">contig%zu\n", i);
      supernode_print(stdout, &db_graph, nodes, orients, len);
      putc('\n', stdout);
    }

    dead_ends += (edges_get_outdegree(db_graph.col_edges[wlk.node], wlk.orient) == 0);
    lengths[i] = len;
    junctions[i] = njunc;
    total_len += len;
    total_junc += njunc;

    density = (double)njunc / len;
    max_density = MAX2(max_density, density);
    max_len = MAX2(max_len, len);
    max_junc = MAX2(max_junc, njunc);
  }

  free(nodes);
  free(orients);

  char total_len_str[100], total_junc_str[100];
  ulong_to_str(total_len, total_len_str);
  ulong_to_str(total_junc, total_junc_str);

  status("\n");
  status("total_len: %s; total_junc: %s (%.2f%% junctions)\n",
         total_len_str, total_junc_str, (100.0*total_junc)/total_len);
  status("dead ends: %zu / %zu\n", dead_ends, num_samples);
  status("mean length: %.2f\n", (double)total_len / num_samples);
  status("mean junctions: %.1f per contig, %.2f%% nodes (1 every %.1f nodes)\n",
          (double)total_junc / num_samples, (100.0 * total_junc) / total_len,
          (double)total_len / total_junc);

  qsort(lengths, num_samples, sizeof(size_t), cmp_size);
  qsort(junctions, num_samples, sizeof(size_t), cmp_size);

  double median_len = MEDIAN(lengths, num_samples);
  double median_junc = MEDIAN(junctions, num_samples);

  status("Median contig length: %.2f\n", median_len);
  status("Median junctions per contig: %.2f\n", median_junc);
  status("Longest contig length: %zu\n", max_len);
  status("Most junctions: %zu\n", max_junc);
  status("Highest junction density: %.2f\n", max_density);
  status("Contig ends which loop %zu [%.2f%%]\n", nloop,
         (100.0 * nloop) / (2.0 * num_samples));

  free(visited);
  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free((void*)db_graph.kmer_paths);
  free(path_store);

  paths_header_dealloc(&pheader);
  graph_walker_dealloc(&wlk);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&file);

  return EXIT_SUCCESS;
}
