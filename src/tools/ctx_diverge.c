#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "path_store.h"
#include "seq_reader.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_walker.h"

static const char usage[] =
"usage: "CMD" diverge [-m <mem>] <in.ctx> <colour> <ref.fa> <out.bubbles.gz>\n"
"  Make bubble calls using a trusted path\n";

typedef struct
{
  uint32_t pos, prev;
} LinkedChromPos;

typedef struct
{
  LinkedChromPos *linkedlist;
  size_t list_len, list_cap;
  uint32_t *kmer_pos;
  GraphWalker wlk;
} DivergeData;


#define getpos(data,node,orient) ((data)->kmer_pos[2*(node)+(orient)])

// Method copied from seq_reader.c
static void load_chrom(const read_t *r, dBGraph *db_graph,
                       int qual_cutoff, int hp_cutoff,
                       Colour colour, SeqLoadingStats *stats,
                       DivergeData *data)
{
  const uint32_t kmer_size = db_graph->kmer_size;
  if(r->seq.end < kmer_size) {
    stats->total_bad_reads++;
    return;
  }

  size_t search_start = 0;
  size_t contig_start, contig_end = 0;

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qual_cutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qual_cutoff, hp_cutoff, &search_start);

    size_t i, contig_len = contig_end - contig_start;

    // printf("contig: %.*s\n", (int)contig_len, r->seq.b+contig_start);

    // Load into graph
    BinaryKmer bkmer, tmp_key;
    hkey_t prev_node, node;
    Orientation prev_or, or;

    bkmer = binary_kmer_from_str(r->seq.b+contig_start, kmer_size);
    tmp_key = db_node_get_key(bkmer, kmer_size);
    prev_node = db_graph_find_or_add_node(db_graph, tmp_key, colour);
    prev_or = db_node_get_orientation(bkmer, tmp_key);

    for(i = contig_start+kmer_size; i < contig_end; i++)
    {
      Nucleotide nuc = binary_nuc_from_char(r->seq.b[i]);
      binary_kmer_left_shift_add(&bkmer, kmer_size, nuc);

      tmp_key = db_node_get_key(bkmer, kmer_size);
      node = db_graph_find_or_add_node(db_graph, tmp_key, colour);
      or = db_node_get_orientation(bkmer, tmp_key);

      db_graph_add_edge(db_graph, colour, prev_node, node, prev_or, or);

      size_t pos = i - kmer_size + 1;
      LinkedChromPos newpos = {.pos = pos, .prev = getpos(data,node,or)};
      memcpy(data->linkedlist + data->list_len, &newpos, sizeof(LinkedChromPos));
      getpos(data,node,or) = data->list_len;
      data->list_len++;

      prev_node = node;
      prev_or = or;
    }

    // Update contig stats
    if(stats->readlen_count_array != NULL) {
      contig_len = MIN2(contig_len, stats->readlen_count_array_size-1);
      stats->readlen_count_array[contig_len]++;
    }
    stats->total_bases_loaded += contig_len;
    stats->kmers_loaded += contig_len + 1 - kmer_size;
  }

  // contig_end == 0 if no contigs from this read
  if(contig_end == 0) stats->total_bad_reads++;
  else stats->total_good_reads++;
}

// Attempt to reconnect with the ref
// Stop at the first node in the ref
static void diverge_call_path(hkey_t node, Orientation orient,
                              GraphWalker *wlk)
{

}

static void diverge_call_node(BinaryKmer bkmer, const dBGraph *db_graph,
                              DivergeData *data)
{
  // DEV: check if node already been used

  GraphWalker *wlk = &data->wlk;

  BinaryKmer bkey = db_node_get_key(bkmer, db_graph->kmer_size);
  hkey_t node = hash_table_find(&db_graph->ht, bkey);
  Orientation orient = db_node_get_orientation(bkmer, bkey);

  // Check for fork in pop and not in ref
  Edges col0edges = db_graph->edges[node] &~ db_node_col_edges(db_graph, 0, node);

  Edges edges;
  Nucleotide nuc;
  Colour colour;

  for(orient = 0; orient < 2; orient++) {
    edges = edges_with_orientation(col0edges, orient);
    for(nuc = 0; nuc < 4; nuc++) {
      if(edges & nuc) {
        for(colour = 1; colour < db_graph->num_of_cols_used; colour++) {
          graph_walker_init(wlk, db_graph, colour, node, orient);

          // DEV: call path
          diverge_call_path(node, orient, wlk);

          graph_walker_finish(wlk);
        }
      }
    }
  }
}

void diverge_call(read_t *r1, read_t *r2,
                  int qoffset1, int qoffset2,
                  const SeqLoadingPrefs *prefs,
                  SeqLoadingStats *stats,
                  void *ptr)
{
  (void)r2; // always NULL for path divergence calling
  (void)qoffset1; (void)qoffset2;

  dBGraph *db_graph = prefs->db_graph;
  const uint32_t kmer_size = db_graph->kmer_size;

  DivergeData *data = (DivergeData*)ptr;

  if(r1->seq.end < kmer_size) return;
  if(r1->seq.end > data->list_cap) {
    data->list_cap = ROUNDUP2POW(r1->seq.end);
    data->linkedlist = realloc2(data->linkedlist,
                               data->list_cap * sizeof(LinkedChromPos));
  }

  memset(data->kmer_pos, 0xff, 2 * db_graph->ht.capacity);

  Colour refcol = 0;
  load_chrom(r1, db_graph, 0, 0, refcol, stats, data);

  READ_TO_BKMERS(r1, kmer_size, 0, 0, stats,
                 diverge_call_node, db_graph, data);

  // DEV: find replacement
  // db_graph_wipe_colour(db_graph, 0);
}

// DEV: load ref into colour 0
//      load pop into colours 1...N
int ctx_diverge(CmdArgs *args)
{
  cmd_accept_options(args, "m");
  cmd_require_options(args, "m", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 4) print_usage(usage, NULL);

  size_t mem_to_use = args->mem_to_use;

  char *input_ctx_path, *input_fa_path, *output_bubble_path;
  uint32_t colour;

  if(!mem_to_integer(argv[1], &mem_to_use) || mem_to_use == 0)
    print_usage(usage, "Invalid memory argument: %s", argv[1]);

  input_ctx_path = argv[2];

  // Probe ctx
  boolean is_binary = false;
  GraphFileHeader gheader = {.capacity = 0};

  if(!graph_file_probe(input_ctx_path, &is_binary, &gheader))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  uint32_t kmer_size = gheader.kmer_size, num_of_cols = gheader.num_of_cols;

  if(!parse_entire_uint(argv[3], &colour) || colour == 0)
    print_usage(usage, "Invalid colour: %s", argv[3]);

  input_fa_path = argv[4];
  if(!test_file_readable(input_fa_path))
    print_usage(usage, "Cannot read trusted reference: %s", input_fa_path);

  output_bubble_path = argv[5];
  if(!test_file_writable(output_bubble_path))
    print_usage(usage, "Cannot write to output file: %s", output_bubble_path);

  // DEV: look for ctp file
  size_t path_mem = 0;

  // Data required:
  // hashtable
  // Edges in graph samples (union)
  // Edges in ref
  // node_in_cols + ref
  // visited fw/rv
  // chrom kmer pos

  // Decide on memory
  size_t hash_kmers, req_num_kmers = gheader.num_of_kmers / IDEAL_OCCUPANCY;
  size_t hash_mem = hash_table_mem(req_num_kmers, &hash_kmers);

  size_t graph_mem = hash_mem +
                     hash_kmers * sizeof(Edges) * 2 + // col_edges
                     hash_kmers * sizeof(uint64_t) + // kmer_paths
                     round_bits_to_bytes(hash_kmers) * num_of_cols + // in col
                     round_bits_to_bytes(hash_kmers) * 2 + // visited fw/rv
                     round_bits_to_bytes(hash_kmers) + // called from node
                     hash_kmers * sizeof(uint32_t) * 2 + // ref coords
                     (0x1<<28) * sizeof(LinkedChromPos);


  // memory to strings
  char graph_mem_str[100], path_mem_str[100], total_mem_str[100];
  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(path_mem, 1, path_mem_str);
  bytes_to_str(graph_mem+path_mem, 1, total_mem_str);

  if(graph_mem > mem_to_use)
    print_usage(usage, "Not enough memory; requires %s", total_mem_str);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, num_of_cols, hash_kmers);

  message("[memory]  graph: %s;  paths: %s\n", graph_mem_str, path_mem_str);

  // Edges
  db_graph.col_edges = calloc2(hash_kmers*2, sizeof(uint8_t));

  // In colour - used is traversal
  size_t words64_per_col = round_bits_to_words64(hash_kmers);
  db_graph.node_in_cols = calloc2(words64_per_col*num_of_cols, sizeof(uint64_t));

  // Paths
  db_graph.kmer_paths = malloc2(hash_kmers * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, hash_kmers * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  path_store_init(&db_graph.pdata, path_store, path_mem, num_of_cols);

  // Allocate memory for calling
  uint32_t *kmer_pos = malloc2(2 * db_graph.ht.capacity * sizeof(uint32_t));
  LinkedChromPos *linkedlist = malloc2((0x1<<28) * sizeof(LinkedChromPos));
  GraphWalker wlk;
  graph_walker_alloc(&wlk);
  DivergeData data = {.linkedlist = linkedlist,
                      .list_len = 0, .list_cap = (0x1<<28),
                      .kmer_pos = kmer_pos, .wlk = wlk};

  // Load graph
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .merge_colours = false,
                           .boolean_covgs = false,
                           .load_seq = false,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = true,
                           .must_exist_in_graph = false,
                           .empty_colours = false,
                           .db_graph = &db_graph};

  graph_load(input_ctx_path, &prefs, stats, NULL);
  hash_table_print_stats(&db_graph.ht);

  read_t r1;
  seq_read_alloc(&r1);
  prefs.load_seq = true;

  seq_parse_se(input_fa_path, &r1, NULL, &prefs, stats, diverge_call, &data);

  seq_read_dealloc(&r1);
  graph_walker_dealloc(&wlk);
  free(kmer_pos);
  free(linkedlist);

  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free((void *)db_graph.kmer_paths);
  free(path_store);

  graph_header_dealloc(&gheader);
  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  message("  Contigs written to: %s\n", output_bubble_path);
  message("Done.\n");
  return EXIT_SUCCESS;
}
