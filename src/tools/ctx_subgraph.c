#include "global.h"

#include "string_buffer.h"
#include "seq_file.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "seq_reader.h"


static const char usage[] =
"usage: "CMD" subgraph [options] <out.ctx> <dist> <in.ctx>[:cols] [in2.ctx ...]\n"
"\n"
"  Loads graphs (in.ctx) and dumps a graph (out.ctx) that contains all kmers within\n"
"  <dist> edges of kmers in <seeds.fa>.  Maintains number of colours / covgs etc.\n"
"\n"
"  Options:\n"
"    -m <mem>          Memory to use\n"
"    -h <kmers>        Hash size\n"
"    --seed <seed.fa>  Read in a seed file\n";

typedef struct
{
  hkey_t *nodes;
  size_t len, capacity;
} dBNodeList;


dBGraph db_graph;
uint64_t *kmer_mask;

static void mark_bkmer(BinaryKmer bkmer, SeqLoadingStats *stats)
{
  #ifdef DEBUG
    char tmp[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph.kmer_size, tmp);
    message("got bkmer %s\n", tmp);
  #endif

  BinaryKmer bkey = db_node_get_key(bkmer, db_graph.kmer_size);
  hkey_t node = hash_table_find(&db_graph.ht, bkey);
  if(node != HASH_NOT_FOUND) bitset_set(kmer_mask, node);
  else stats->unique_kmers++;
}

void mark_reads(read_t *r1, read_t *r2,
                int qoffset1, int qoffset2,
                const SeqLoadingPrefs *prefs, SeqLoadingStats *stats, void *ptr)
{
  (void)qoffset1;
  (void)qoffset2;
  (void)prefs;
  (void)ptr;

  READ_TO_BKMERS(r1, db_graph.kmer_size, 0, 0, stats, mark_bkmer, stats);
  if(r2 != NULL) {
    READ_TO_BKMERS(r2, db_graph.kmer_size, 0, 0, stats, mark_bkmer, stats);
  }
}

static void store_node_neighbours(const hkey_t node, dBNodeList *list)
{
  // Get neighbours
  Edges edges = db_graph.col_edges[node];
  int num_next, i;
  hkey_t next_nodes[8];
  Orientation next_orients[8];
  Nucleotide next_bases[8];

  BinaryKmer bkmer = db_node_bkmer(&db_graph, node);

  // Get neighbours in forward dir
  num_next  = db_graph_next_nodes(&db_graph, bkmer, FORWARD, edges,
                                  next_nodes, next_orients, next_bases);

  // Get neighbours in reverse dir
  num_next += db_graph_next_nodes(&db_graph, bkmer, REVERSE, edges,
                                  next_nodes+num_next, next_orients+num_next,
                                  next_bases+num_next);

  // if not flagged add to list
  for(i = 0; i < num_next; i++) {
    if(!bitset_has(kmer_mask, next_nodes[i])) {
      bitset_set(kmer_mask, next_nodes[i]);
      // if list full, exit
      if(list->len == list->capacity) die("Please increase <mem> size");
      list->nodes[list->len++] = next_nodes[i];
    }
  }
}

static void store_bkmer_neighbours(BinaryKmer bkmer, dBNodeList *list,
                                   SeqLoadingStats *stats)
{
  BinaryKmer bkey = db_node_get_key(bkmer, db_graph.kmer_size);
  hkey_t node = hash_table_find(&db_graph.ht, bkey);
  if(node != HASH_NOT_FOUND) store_node_neighbours(node, list);
  else stats->unique_kmers++;
}

void store_nodes(read_t *r1, read_t *r2,
                 int qoffset1, int qoffset2,
                 const SeqLoadingPrefs *prefs, SeqLoadingStats *stats, void *ptr)
{
  (void)qoffset1;
  (void)qoffset2;
  (void)prefs;

  dBNodeList *list = (dBNodeList*)ptr;

  READ_TO_BKMERS(r1, db_graph.kmer_size, 0, 0, stats,
                 store_bkmer_neighbours, list, stats);

  if(r2 != NULL) {
    READ_TO_BKMERS(r2, db_graph.kmer_size, 0, 0, stats,
                   store_bkmer_neighbours, list, stats);
  }
}

int ctx_subgraph(CmdArgs *args)
{
  cmd_accept_options(args, "mh");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 4) print_usage(usage, NULL);

  char *seed_files[argc];
  size_t num_seed_files = 0;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcasecmp(argv[argi],"--seed") == 0)
    {
      if(argi+1 == argc)
        print_usage(usage, "--seed <seed.fa> requires and argument");
      seed_files[num_seed_files++] = argv[argi+1];
      if(!test_file_readable(argv[argi+1]))
        die("Cannot read seed file: %s", argv[argi+1]);
      argi++;
    }
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  char *out_path = argv[argi], *diststr = argv[argi+1];
  uint32_t dist;

  if(!parse_entire_uint(diststr, &dist))
    print_usage(usage, "Invalid <dist> value, must be int >= 0: %s", diststr);

  int num_binaries_int = argc - 2*num_seed_files - 2;
  if(num_binaries_int <= 0)
    print_usage(usage, "Please specify input graph files (.ctx)");

  size_t i, num_binaries = num_binaries_int;
  char **binary_paths = argv + 2*num_seed_files + 2;

  //
  // Probe binaries to get kmer_size
  //
  boolean is_binary = false;
  uint32_t kmer_size;
  uint32_t ctx_num_of_cols[num_binaries], ctx_max_cols[num_binaries];
  uint64_t max_num_kmers = 0;
  GraphFileHeader gheader = {.capacity = 0};

  for(i = 0; i < num_binaries; i++)
  {
    if(!graph_file_probe(binary_paths[i], &is_binary, &gheader))
      print_usage(usage, "Cannot read input binary file: %s", binary_paths[i]);
    else if(!is_binary)
      print_usage(usage, "Input binary file isn't valid: %s", binary_paths[i]);

    if(i == 0) kmer_size = gheader.kmer_size;
    else if(kmer_size != gheader.kmer_size) {
      die("Graph kmer-sizes do not match [%u vs %u; %s; %s]\n",
          kmer_size, gheader.kmer_size, binary_paths[i-1], binary_paths[i]);
    }

    ctx_num_of_cols[i] = gheader.num_of_cols;
    ctx_max_cols[i] = gheader.max_col;
    max_num_kmers = MAX2(gheader.num_of_kmers, max_num_kmers);
  }

  //
  // Calculate memory use
  //
  size_t kmers_in_hash, ideal_capacity = max_num_kmers / IDEAL_OCCUPANCY;
  size_t req_num_kmers = args->num_kmers_set ? args->num_kmers : ideal_capacity;
  size_t hash_mem = hash_table_mem(req_num_kmers, &kmers_in_hash);
  size_t kmer_mask_mem = round_bits_to_words64(kmers_in_hash) * sizeof(uint64_t);
  size_t fringe_mem = args->mem_to_use - hash_mem - kmer_mask_mem;
  size_t num_of_fringe_nodes = fringe_mem / (sizeof(hkey_t) * 2);
  size_t search_mem = num_of_fringe_nodes * (sizeof(hkey_t) * 2);

  char num_kmers_str[100];
  ulong_to_str(max_num_kmers, num_kmers_str);

  char kmers_in_hash_str[100], search_mem_str[100];
  bytes_to_str(kmers_in_hash, 1, kmers_in_hash_str);
  bytes_to_str(search_mem, 1, search_mem_str);

  message("[memory]  graph: %s; search: %s\n", kmers_in_hash_str, search_mem_str);

  if(hash_mem + kmer_mask_mem > args->mem_to_use) {
    char mem_str[100];
    bytes_to_str(hash_mem + kmer_mask_mem, 1, mem_str);
    print_usage(usage, "Require more memory (-m <mem>) [suggested > %s]", mem_str);
  }
  else if(kmers_in_hash < max_num_kmers) {
    print_usage(usage, "Not enough kmers in the hash, require: %s "
                       "(set bigger -h <kmers> or -m <mem>)", num_kmers_str);
  }
  else if(kmers_in_hash < max_num_kmers / WARN_OCCUPANCY)
    warn("Low memory for binary size (require: %s)", num_kmers_str);

  if(num_of_fringe_nodes < 100)
    die("Not enough memory for the graph search (set -m <mem> higher)");

  if(!test_file_writable(out_path))
    die("Cannot write to output file: %s", out_path);

  // Create db_graph with one colour
  db_graph_alloc(&db_graph, kmer_size, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity, sizeof(Covg));

  size_t num_words64 = round_bits_to_words64(db_graph.ht.capacity);
  kmer_mask = calloc2(num_words64, sizeof(uint64_t));

  // Store edge nodes here
  dBNodeList list0, list1, listtmp;
  list0.nodes = malloc2(sizeof(hkey_t) * num_of_fringe_nodes);
  list1.nodes = malloc2(sizeof(hkey_t) * num_of_fringe_nodes);
  list0.capacity = list1.capacity = num_of_fringe_nodes;
  list0.len = list1.len = 0;

  //
  // Load binaries
  //
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                           .boolean_covgs = false,
                           .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_graph = false,
                           .empty_colours = true,
                           .db_graph = &db_graph};

  for(i = 0; i < num_binaries; i++)
    graph_load(binary_paths[i], &prefs, stats, NULL);

  size_t num_of_binary_kmers = stats->kmers_loaded;

  // Load sequence and mark in first pass
  read_t r1, r2;
  seq_read_alloc(&r1);
  seq_read_alloc(&r2);
  for(i = 0; i < num_seed_files; i++)
    seq_parse_se(seed_files[i], &r1, &r2, &prefs, stats, mark_reads, NULL);
  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  size_t num_of_seed_kmers = stats->kmers_loaded - num_of_binary_kmers;

  message("Read in %zu seed kmers\n", num_of_seed_kmers);

  hash_table_print_stats(&db_graph.ht);

  if(dist > 0)
  {
    // Get edge nodes
    for(i = 0; i < num_seed_files; i++)
      seq_parse_se(seed_files[i], &r1, &r2, &prefs, stats, store_nodes, &list0);

    size_t i, d;
    for(d = 1; d < dist; d++)
    {
      for(i = 0; i < list0.len; i++) {
        store_node_neighbours(list0.nodes[i], &list1);
      }
      list0.len = 0;
      SWAP(list0, list1, listtmp);
    }
  }

  // free
  free(list0.nodes);
  free(list1.nodes);

  message("Pruning untouched nodes...\n");

  // Remove nodes that were not flagged
  db_graph_prune_nodes_lacking_flag(&db_graph, kmer_mask);
  hash_table_print_stats(&db_graph.ht);

  // Dump nodes that were flagged
  if(num_binaries == 1 && ctx_num_of_cols[0] == 1)
  {
    // We have all the info to dump now
    graph_file_save(out_path, &db_graph, CTX_GRAPH_FILEFORMAT, NULL, 0, 1);
  }
  else
  {
    graph_files_merge(out_path, binary_paths, num_binaries,
                      ctx_num_of_cols, ctx_max_cols,
                      false, false, true, &db_graph);
  }

  message("Done.\n");

  free(kmer_mask);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);

  seq_loading_stats_free(stats);
  graph_header_dealloc(&gheader);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
