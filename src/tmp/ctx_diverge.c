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
  dBGraph *db_graph;
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
      Nucleotide nuc = dna_char_to_nuc(r->seq.b[i]);
      bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);

      tmp_key = db_node_get_key(bkmer, kmer_size);
      node = db_graph_find_or_add_node(db_graph, tmp_key, colour);
      or = db_node_get_orientation(bkmer, tmp_key);

      db_graph_add_edge(db_graph, 0, prev_node, node, prev_or, or);

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
  hkey_t hkey = hash_table_find(&db_graph->ht, bkey);
  Orientation orient = db_node_get_orientation(bkmer, bkey);

  // Check for fork in pop and not in ref
  Edges col0edges = db_node_get_edges(db_graph, 1, hkey) &~
                    db_node_get_edges(db_graph, 0, hkey);

  Edges edges;
  Nucleotide nuc;
  Colour colour;

  for(orient = 0; orient < 2; orient++) {
    edges = edges_with_orientation(col0edges, orient);
    for(nuc = 0; nuc < 4; nuc++) {
      if(edges & nuc) {
        for(colour = 1; colour < db_graph->num_of_cols_used; colour++) {
          dBNode node = {.key = hkey, .orient = orient};
          graph_walker_init(wlk, db_graph, colour, colour, node);

          // DEV: call path
          diverge_call_path(hkey, orient, wlk);

          graph_walker_finish(wlk);
        }
      }
    }
  }
}

void diverge_call(read_t *r1, read_t *r2,
                  uint8_t qoffset1, uint8_t qoffset2,
                  void *ptr)
{
  (void)r2; // always NULL for path divergence calling
  (void)qoffset1; (void)qoffset2;

  DivergeData *data = (DivergeData*)ptr;
  dBGraph *db_graph = data->db_graph;
  const uint32_t kmer_size = db_graph->kmer_size;

  if(r1->seq.end < kmer_size) return;
  if(r1->seq.end > data->list_cap) {
    data->list_cap = roundup2pow(r1->seq.end);
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
  cmd_accept_options(args, "m", usage);
  cmd_require_options(args, "m", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 4) print_usage(usage, NULL);

  size_t mem_to_use = args->mem_to_use;

  char *input_ctx_path, *output_bubble_path;
  uint32_t colour;
  seq_file_t *input_fa_file;

  if(!mem_to_integer(argv[1], &mem_to_use) || mem_to_use == 0)
    print_usage(usage, "Invalid memory argument: %s", argv[1]);

  input_ctx_path = argv[2];

  // Probe ctx
  boolean is_binary = false;
  GraphFileHeader gheader = INIT_GRAPH_FILE_HDR;

  if(!graph_file_probe(input_ctx_path, &is_binary, &gheader))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  uint32_t kmer_size = gheader.kmer_size, num_of_cols = gheader.num_of_cols;

  if(!parse_entire_uint(argv[3], &colour) || colour == 0)
    print_usage(usage, "Invalid colour: %s", argv[3]);

  if((input_fa_file = seq_open(argv[4])) == NULL)
    print_usage(usage, "Cannot read trusted reference: %s", argv[4]);

  output_bubble_path = argv[5];
  if(!futil_is_file_writable(output_bubble_path))
    print_usage(usage, "Cannot write to output file: %s", output_bubble_path);

  // DEV: look for ctp file
  // size_t path_mem = 0;

  // Data required:
  // hashtable
  // Edges in graph samples (union)
  // Edges in ref
  // node_in_cols + ref
  // visited fw/rv
  // chrom kmer pos

  // load ref into colour 0, population into colour 1

  // DEV: load paths
  size_t path_mem = 0;

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;
  size_t graph_mem, chrom_pos_mem, chrom_pos_num, total_mem;
  char chrom_pos_mem_str[100], chrom_pos_num_str[100];
  char path_mem_str[100];

  // 2 colours: col_edges, kmer_paths, in_col, visited fw/rv, called from node,
  //            ref coords
  bits_per_kmer = sizeof(Edges)*2*8 + sizeof(uint64_t)*8 + num_of_cols + 2 + 1 +
                  sizeof(uint32_t)*2*8;

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        gheader.num_of_kmers, false, &graph_mem);

  chrom_pos_num = (0x1<<28);
  chrom_pos_mem = chrom_pos_num * sizeof(LinkedChromPos);
  bytes_to_str(chrom_pos_mem, 1, chrom_pos_mem_str);
  ulong_to_str(chrom_pos_num, chrom_pos_num_str);

  total_mem = graph_mem + chrom_pos_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  status("[memory] chrom positions: %s (%s) paths: %s\n",
         chrom_pos_num_str, chrom_pos_mem_str, path_mem_str);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, num_of_cols, 2, kmers_in_hash);

  // status("[memory]  graph: %s;  paths: %s\n", graph_mem_str, path_mem_str);

  // Edges
  db_graph.col_edges = calloc2(kmers_in_hash*2, sizeof(uint8_t));

  // In colour - used is traversal
  size_t bytes_per_col = roundup_bits2bytes(kmers_in_hash);
  db_graph.node_in_cols = calloc2(bytes_per_col*num_of_cols, sizeof(uint8_t));

  // Paths
  db_graph.kmer_paths = malloc2(kmers_in_hash * sizeof(PathIndex));
  memset(db_graph.kmer_paths, 0xff, kmers_in_hash * sizeof(PathIndex));

  path_store_alloc(&db_graph.pdata, path_mem, 0, num_of_cols);

  // Allocate memory for calling
  uint32_t *kmer_pos = malloc2(2 * db_graph.ht.capacity * sizeof(uint32_t));
  LinkedChromPos *linkedlist = malloc2((0x1<<28) * sizeof(LinkedChromPos));
  GraphWalker wlk;
  graph_walker_alloc(&wlk);
  DivergeData data = {.db_graph = db_graph, .linkedlist = linkedlist,
                      .list_len = 0, .list_cap = (0x1<<28),
                      .kmer_pos = kmer_pos, .wlk = wlk};

  // Loading
  SeqLoadingStats *stats = seq_loading_stats_create(0);

  // Graph loading prefs
  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = false};

  graph_load(input_ctx_path, gprefs, stats, NULL);
  hash_table_print_stats(&db_graph.ht);

  read_t r1;
  if(seq_read_alloc(&r1) == NULL) die("Out of memory");

  seq_parse_se_sf(input_fa_file, 0, &r1, NULL, diverge_call, &data);

  seq_read_dealloc(&r1);
  graph_walker_dealloc(&wlk);
  free(kmer_pos);
  free(linkedlist);

  free(db_graph.col_edges);
  free(db_graph.node_in_cols);
  free(db_graph.kmer_paths);

  graph_header_dealloc(&gheader);
  seq_loading_stats_free(stats);
  path_store_dealloc(&db_graph.pdata);
  db_graph_dealloc(&db_graph);

  status("Contigs written to: %s\n", output_bubble_path);

  return EXIT_SUCCESS;
}
