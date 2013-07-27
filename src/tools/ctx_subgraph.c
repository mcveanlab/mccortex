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
#include "binary_format.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" <in.ctx> <seeds.falist> <dist> <out.ctx>\n"
"  Loads <in.ctx> and dumps a binary <out.ctx> that contains all kmers within\n"
"  <dist> edges of kmers in <filelist>.  Maintains number of colours / covgs etc.\n"
"  <mem> specifies how much memory to use to store the list of edge kmers.  Will\n"
"  fail if this is too small.  We suggest 1GB if available on your machine.\n";

typedef struct
{
  hkey_t *nodes;
  size_t len, capacity;
} dBNodeList;


dBGraph db_graph;
size_t num_kmers_read = 0;

uint64_t *kmer_mask;

static void mark_bkmer(const BinaryKmer bkmer, SeqLoadingStats *stats)
{
  #ifdef DEBUG
    char tmp[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph.kmer_size, tmp);
    message("got bkmer %s\n", tmp);
  #endif

  BinaryKmer tmpkey;
  db_node_get_key(bkmer, db_graph.kmer_size, tmpkey);
  hkey_t node = hash_table_find(&db_graph.ht, tmpkey);
  if(node != HASH_NOT_FOUND) bitset_set(kmer_mask, node);
  else stats->unique_kmers++;
}

void mark_reads(read_t *r1, read_t *r2,
                int qoffset1, int qoffset2,
                SeqLoadingPrefs *prefs, SeqLoadingStats *stats, void *ptr)
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
  Edges edges = db_graph.edges[node];
  int num_next, i;
  hkey_t next_nodes[8];
  BinaryKmer next_bkmers[8];

  ConstBinaryKmerPtr bkmer = db_node_bkmer(&db_graph, node);

  // Get neighbours in forward dir
  num_next  = db_graph_next_nodes(&db_graph, bkmer, edges & 0xf,
                                  next_nodes, next_bkmers);

  BinaryKmer revbkmer;
  binary_kmer_reverse_complement(bkmer, db_graph.kmer_size, revbkmer);

  // Get neighbours in reverse dir
  num_next += db_graph_next_nodes(&db_graph, revbkmer, edges>>4,
                                  next_nodes+num_next,
                                  next_bkmers+num_next);

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

static void store_bkmer_neighbours(const BinaryKmer bkmer, dBNodeList *list,
                                   SeqLoadingStats *stats)
{
  BinaryKmer tmpkey;
  db_node_get_key(bkmer, db_graph.kmer_size, tmpkey);
  hkey_t node = hash_table_find(&db_graph.ht, tmpkey);
  if(node != HASH_NOT_FOUND) store_node_neighbours(node, list);
  else stats->unique_kmers++;
}

void store_nodes(read_t *r1, read_t *r2,
                 int qoffset1, int qoffset2,
                 SeqLoadingPrefs *prefs, SeqLoadingStats *stats, void *ptr)
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

static void filter_subgraph(const char *input_ctx_path,
                            const char *input_filelist, int dist,
                            size_t num_of_fringe_nodes,
                            const char *out_path)
{
  // Store edge nodes here
  dBNodeList list0, list1, listtmp;
  list0.nodes = malloc(sizeof(hkey_t) * num_of_fringe_nodes);
  list1.nodes = malloc(sizeof(hkey_t) * num_of_fringe_nodes);
  list0.capacity = list1.capacity = num_of_fringe_nodes;
  list0.len = list1.len = 0;

  if(list0.nodes == NULL || list1.nodes == NULL) die("Out of memory");

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = true,
                           .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_colour = -1,
                           .empty_colours = true,
                           .update_ginfo = false,
                           .db_graph = &db_graph};

  // Load binary
  binary_load(input_ctx_path, &db_graph, &prefs, stats, NULL);

  size_t num_of_binary_kmers = stats->kmers_loaded;

  // Load sequence and mark in first pass
  parse_filelists(input_filelist, NULL, READ_FALIST, &prefs, stats,
                  mark_reads, NULL);

  size_t num_of_seed_kmers = stats->kmers_loaded - num_of_binary_kmers;

  if(dist > 0)
  {
    // Get edge nodes
    parse_filelists(input_filelist, NULL, READ_FALIST, &prefs, stats,
                    store_nodes, &list0);  

    int d;
    size_t i;

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

  // Remove nodes that were not flagged
  db_graph_prune_nodes_lacking_flag(&db_graph, kmer_mask);

  // Dump nodes that were flagged
  size_t nodes_dumped = db_graph_filter_file(&db_graph, input_ctx_path, out_path);

  message("Read in %zu seed kmers\n", num_of_seed_kmers);
  message("Dumped %zu kmers\n", nodes_dumped);
  message("Done.\n");

  seq_loading_stats_free(stats);
}

int ctx_subgraph(CmdArgs *args)
{
  cmd_accept_options(args, "m");
  cmd_require_options(args, "m");
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 4) print_usage(usage, NULL);

  uint64_t mem_to_use = args->mem_to_use;

  char *input_ctx_path, *input_filelist, *out_path;
  uint32_t dist;

  input_ctx_path = argv[2];

  input_filelist = argv[3];
  if(!test_file_readable(input_filelist))
    print_usage(usage, "Cannot read filelist: %s", input_filelist);

  if(!parse_entire_uint(argv[4], &dist))
    print_usage(usage, "Invalid <dist> value -- must be int >= 0: %s", argv[4]);

  out_path = argv[5];

  // Probe binary to get kmer_size
  boolean is_binary = false;
  uint32_t kmer_size, num_of_cols, max_col;
  uint64_t num_kmers;

  if(!binary_probe(input_ctx_path, &is_binary, &kmer_size, &num_of_cols, &max_col, &num_kmers))
    print_usage(usage, "Cannot read input binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  size_t num_of_hash_kmers;
  size_t req_num_kmers = num_kmers*(1.0/IDEAL_OCCUPANCY);
  size_t hash_mem = hash_table_mem(req_num_kmers, &num_of_hash_kmers);
  size_t fringe_mem = mem_to_use - hash_mem;
  size_t num_of_fringe_nodes = fringe_mem / (sizeof(hkey_t) * 2);

  if(hash_mem >= mem_to_use || num_of_fringe_nodes < 100)
    die("Not enough memory for the graph");

  if(!test_file_writable(out_path))
    die("Cannot write to output file: %s", out_path);

  message("Using kmer size: %u\n", kmer_size);
  message("Using %zu bytes for graph search\n",
          num_of_fringe_nodes * (sizeof(hkey_t) * 2));

  // Create db_graph
  db_graph_alloc(&db_graph, kmer_size, 1, num_of_hash_kmers);
  db_graph.edges = calloc(db_graph.ht.capacity, sizeof(Edges));

  size_t num_words64 = round_bits_to_words64(db_graph.ht.capacity);
  kmer_mask = calloc(num_words64, sizeof(uint64_t));

  filter_subgraph(input_ctx_path, input_filelist, dist,
                  num_of_fringe_nodes, out_path);

  free(kmer_mask);
  free(db_graph.edges);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
