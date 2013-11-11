#include "global.h"

#include "string_buffer.h"
#include "seq_file.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "graph_info.h"
#include "db_node.h"
#include "graph_format.h"
#include "seq_reader.h"
#include "prune_nodes.h"

// DEV: add --usecols <c> option
static const char usage[] =
"usage: "CMD" subgraph [options] <out.ctx> <dist> <in.ctx>[:cols] [in2.ctx ...]\n"
"\n"
"  Loads graphs (in.ctx) and dumps a graph (out.ctx) that contains all kmers within\n"
"  <dist> edges of kmers in <seeds.fa>.  Maintains number of colours / covgs etc.\n"
"  Loads seed files twice: 1) get seed; 2) extend;  This lowers memory requirement\n"
"  for large (seed) graphs but means seed files cannot be pipes / sockets.\n"
"\n"
"  Options:\n"
"    -m <mem>          Memory to use  <required>\n"
"    -h <kmers>        Hash size\n"
"    --seed <seed.fa>  Read in a seed file\n"
"    --invert          Dump kmers not in subgraph\n"
"    --usecols <n>     Number of samples in memory at once (speedup)\n";

typedef struct
{
  hkey_t *nodes;
  size_t len, capacity;
} dBNodeList;


dBGraph db_graph;
uint64_t *kmer_mask;

static void mark_bkmer(BinaryKmer bkmer, SeqLoadingStats *stats)
{
  #ifdef CTXVERBOSE
    char tmp[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph.kmer_size, tmp);
    status("got bkmer %s\n", tmp);
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
  Edges edges = db_node_edges_union(&db_graph, node);
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
  boolean invert = false;
  size_t usencols = 1;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcasecmp(argv[argi], "--seed") == 0)
    {
      if(argi+1 == argc)
        print_usage(usage, "--seed <seed.fa> requires an argument");
      seed_files[num_seed_files] = argv[argi+1];
      if(!test_file_readable(seed_files[num_seed_files]))
        die("Cannot read --seed file: %s", argv[argi+1]);
      argi++; num_seed_files++;
    }
    else if(strcasecmp(argv[argi], "--invert") == 0) invert = true;
    else if(strcasecmp(argv[argi], "--usecols") == 0)
    {
      if(argi+1 == argc || !parse_entire_size(argv[argi+1], &usencols) || !usencols)
        print_usage(usage, "--usecols <C> requires a +ve integer argument");
      argi++;
    }
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  char *out_path = argv[argi], *diststr = argv[argi+1];
  uint32_t dist;

  if(!parse_entire_uint(diststr, &dist))
    print_usage(usage, "Invalid <dist> value, must be int >= 0: %s", diststr);

  int num_files_int = argc - 2*num_seed_files - 2;
  if(num_files_int <= 0)
    print_usage(usage, "Please specify input graph files (.ctx)");

  size_t i, j, col, num_files = num_files_int, total_cols = 0;
  char **paths = argv + 2*num_seed_files + 2;

  // Open graph files
  uint64_t max_num_kmers = 0;
  GraphFileReader files[num_files];

  for(i = 0; i < num_files; i++)
  {
    files[i] = INIT_GRAPH_READER;
    graph_file_open(&files[i], paths[i], true);

    if(files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      die("Graph kmer-sizes do not match [%u vs %u; %s; %s]\n",
          files[0].hdr.kmer_size, files[i].hdr.kmer_size,
          files[0].path.buff, files[i].path.buff);
    }

    size_t offset = total_cols;
    total_cols += files[i].intocol + graph_file_outncols(&files[i]);
    files[i].intocol += offset;

    max_num_kmers = MAX2(files[i].hdr.num_of_kmers, max_num_kmers);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t num_of_fringe_nodes, fringe_mem;
  char graph_mem_str[100], num_fringe_nodes_str[100], fringe_mem_str[100];

  bits_per_kmer = ((sizeof(Edges) + sizeof(Covg))*usencols*8 + 1);
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, max_num_kmers, false);

  graph_mem = hash_table_mem(kmers_in_hash,false,NULL) +
              (kmers_in_hash*bits_per_kmer)/8;
  bytes_to_str(graph_mem, 1, graph_mem_str);

  if(graph_mem >= args->mem_to_use)
    die("Not enough memory for graph (requires %s)", graph_mem_str);

  // Fringe nodes
  fringe_mem = args->mem_to_use - graph_mem;
  num_of_fringe_nodes = fringe_mem / (sizeof(hkey_t) * 2);
  bytes_to_str(fringe_mem, 1, fringe_mem_str);
  ulong_to_str(num_of_fringe_nodes, num_fringe_nodes_str);
  status("[memory] fringe nodes: %s (%s)\n", num_fringe_nodes_str, fringe_mem_str);

  if(num_of_fringe_nodes < 100)
    die("Not enough memory for the graph search (set -m <mem> higher)");

  if(!test_file_writable(out_path))
    die("Cannot write to output file: %s", out_path);

  // Create db_graph with one colour
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, usencols, usencols, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity*usencols, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity*usencols, sizeof(Covg));

  size_t num_words64 = round_bits_to_words64(db_graph.ht.capacity*usencols);
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
  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           // binaries
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .must_exist_in_edges = NULL,
                           .empty_colours = true,
                           // seq
                           .into_colour = 0,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false};

  StrBuf intersect_gname;
  strbuf_alloc(&intersect_gname, 1024);

  size_t tmpcol; boolean tmpflatten;
  for(i = 0; i < num_files; i++) {
    tmpcol = files[i].intocol;
    tmpflatten = files[i].flatten;

    if(total_cols > db_graph.num_of_cols) {
      files[i].intocol = 0;
      files[i].flatten = true;
    }

    graph_load(&files[i], &prefs, stats);
    files[i].intocol = tmpcol;
    files[i].flatten = tmpflatten;

    for(j = 0; j < files[i].ncols; j++) {
      col = files[i].cols[j];
      graph_info_make_intersect(&files[i].hdr.ginfo[col], &intersect_gname);
    }
  }

  char subgraphstr[] = "subgraph:{";
  strbuf_insert(&intersect_gname, 0, subgraphstr, strlen(subgraphstr));
  strbuf_append_char(&intersect_gname, '}');

  size_t num_of_binary_kmers = stats->kmers_loaded;

  // Load sequence and mark in first pass
  read_t r1, r2;
  seq_read_alloc(&r1);
  seq_read_alloc(&r2);
  for(i = 0; i < num_seed_files; i++)
    seq_parse_se(seed_files[i], &r1, &r2, &prefs, stats, mark_reads, NULL);

  size_t num_of_seed_kmers = stats->kmers_loaded - num_of_binary_kmers;

  status("Read in %zu seed kmers\n", num_of_seed_kmers);

  hash_table_print_stats(&db_graph.ht);

  if(dist > 0)
  {
    char diststr[100];
    ulong_to_str(dist, diststr);
    status("Extending subgraph by %s\n", diststr);

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

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  free(list0.nodes);
  free(list1.nodes);

  status("Pruning untouched nodes...\n");

  if(invert) {
    for(i = 0; i < num_words64; i++)
      kmer_mask[i] = ~kmer_mask[i];
  }

  // Remove nodes that were not flagged
  prune_nodes_lacking_flag(&db_graph, kmer_mask);
  free(kmer_mask);

  // Dump nodes that were flagged
  Edges *intersect_edges = NULL;
  boolean kmers_loaded = true;
  boolean colours_loaded = (total_cols <= db_graph.num_of_cols);

  if(!colours_loaded)
  {
    // Need to reload graph colours - therefore construct edge intersection set
    intersect_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
    for(i = 0; i < db_graph.ht.capacity; i++)
      intersect_edges[i] = db_node_edges_union(&db_graph, i);
  }

  graph_files_merge_mkhdr(out_path, files, num_files,
                          kmers_loaded, colours_loaded,
                          intersect_edges, intersect_gname.buff,
                          &db_graph);

  if(intersect_edges != NULL) free(intersect_edges);
  free(db_graph.col_edges);
  free(db_graph.col_covgs);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);

  strbuf_dealloc(&intersect_gname);
  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
