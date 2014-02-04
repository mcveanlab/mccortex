#include "global.h"

#include "string_buffer.h"
#include "seq_file.h"

#include "tools.h"
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

const char subgraph_usage[] =
"usage: "CMD" subgraph [options] <out.ctx> <dist> <in.ctx>[:cols] [in2.ctx ...]\n"
"\n"
"  Loads graphs (in.ctx) and dumps a graph (out.ctx) that contains all kmers within\n"
"  <dist> edges of kmers in <seeds.fa>.  Maintains number of colours / covgs etc.\n"
"  Loads seed files twice: 1) get seed; 2) extend;  This lowers memory requirement\n"
"  for large (seed) graphs but means seed files cannot be pipes / sockets.\n"
"\n"
"  Options:\n"
"    -m <mem>          Memory to use  <required>\n"
"    -n <kmers>        Hash size\n"
"    --seq <seed.fa>   Read in a seed file\n"
"    --invert          Dump kmers not in subgraph\n"
"    --ncols <n>       Number of samples in memory at once (speedup)\n";

typedef struct
{
  const dBGraph *const db_graph;
  uint64_t *const kmer_mask;
  dBNodeBuffer nbuf;
  LoadingStats stats;
} EdgeNodes;

static void edge_nodes_alloc(EdgeNodes *enodes, const dBGraph *graph,
                             uint64_t *kmer_mask, size_t capacity)
{
  EdgeNodes tmp = {.db_graph = graph, .kmer_mask = kmer_mask};
  memcpy(enodes, &tmp, sizeof(EdgeNodes));
  db_node_buf_alloc(&enodes->nbuf, capacity);
  loading_stats_init(&enodes->stats);
}

static void edge_nodes_dealloc(EdgeNodes *enodes)
{
  db_node_buf_dealloc(&enodes->nbuf);
}

static void mark_bkmer(BinaryKmer bkmer, EdgeNodes *enodes)
{
  const dBGraph *db_graph = enodes->db_graph;
  dBNode node = db_graph_find(db_graph, bkmer);

  #ifdef CTXVERBOSE
    char tmp[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph->kmer_size, tmp);
    status("got bkmer %s\n", tmp);
  #endif

  if(node.key != HASH_NOT_FOUND) {
    if(!bitset_get(enodes->kmer_mask, node.key) && enodes->nbuf.capacity > 0 &&
       !db_node_buf_attempt_add(&enodes->nbuf, node)) {
      die("Please increase <mem> size");
    }
    bitset_set(enodes->kmer_mask, node.key);
  }
}

void store_read_nodes(read_t *r1, read_t *r2,
                      uint8_t qoffset1, uint8_t qoffset2, void *ptr)
{
  (void)qoffset1; (void)qoffset2;
  EdgeNodes *enodes = (EdgeNodes*)ptr;
  const dBGraph *db_graph = enodes->db_graph;

  READ_TO_BKMERS(r1, db_graph->kmer_size, 0, 0, &enodes->stats, mark_bkmer, enodes);
  if(r2 != NULL) {
    READ_TO_BKMERS(r2, db_graph->kmer_size, 0, 0, &enodes->stats, mark_bkmer, enodes);
  }
}

static void store_node_neighbours(const hkey_t hkey, EdgeNodes *enodes)
{
  // Get neighbours
  const dBGraph *db_graph = enodes->db_graph;
  BinaryKmer bkmer = db_node_get_bkmer(db_graph, hkey);
  Edges edges = db_node_get_edges_union(db_graph, hkey);
  size_t num_next, i;
  dBNode next_nodes[8];
  Nucleotide next_bases[8];

  // Get neighbours in forward dir
  num_next  = db_graph_next_nodes(db_graph, bkmer, FORWARD, edges,
                                  next_nodes, next_bases);

  // Get neighbours in reverse dir
  num_next += db_graph_next_nodes(db_graph, bkmer, REVERSE, edges,
                                  next_nodes+num_next, next_bases+num_next);

  // if not flagged add to list
  for(i = 0; i < num_next; i++) {
    if(!bitset_get(enodes->kmer_mask, next_nodes[i].key) &&
       !db_node_buf_attempt_add(&enodes->nbuf, next_nodes[i])) {
      die("Please increase <mem> size");
    }
    bitset_set(enodes->kmer_mask, next_nodes[i].key);
  }
}


int ctx_subgraph(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Already checked that we have at least 4 args

  char *seed_files[argc];
  size_t num_seed_files = 0;
  boolean invert = false;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(!strcasecmp(argv[argi], "--seq") | !strcasecmp(argv[argi], "--seed"))
    {
      if(argi+1 == argc)
        cmd_print_usage("%s <seed.fa> requires an argument", argv[argi]);
      seed_files[num_seed_files] = argv[argi+1];
      if(!futil_is_file_readable(seed_files[num_seed_files]))
        die("Cannot read %s file: %s", argv[argi], argv[argi+1]);
      argi++; num_seed_files++;
    }
    else if(strcasecmp(argv[argi], "--invert") == 0) invert = true;
    else cmd_print_usage("Unknown option: %s", argv[argi]);
  }

  const char *out_path = argv[argi], *diststr = argv[argi+1];
  uint32_t dist;

  if(!parse_entire_uint(diststr, &dist))
    cmd_print_usage("Invalid <dist> value, must be int >= 0: %s", diststr);

  int num_files_int = argc - 2*(int)num_seed_files - 2;
  if(num_files_int <= 0)
    cmd_print_usage("Please specify input graph files (.ctx)");

  size_t i, j, col, num_files = (size_t)num_files_int, total_cols = 0;
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
          files[0].fltr.file_path.buff, files[i].fltr.file_path.buff);
    }

    size_t offset = total_cols;
    total_cols += graph_file_usedcols(&files[i]);
    file_filter_update_intocol(&files[i].fltr, files[i].fltr.intocol + offset);

    max_num_kmers = MAX2(files[i].hdr.num_of_kmers, max_num_kmers);
  }

  //
  // Decide on memory
  //
  const size_t use_ncols = MIN2(args->use_ncols, total_cols);
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t num_of_fringe_nodes = 0, fringe_mem;
  char graph_mem_str[100], num_fringe_nodes_str[100], fringe_mem_str[100];

  bits_per_kmer = ((sizeof(Edges) + sizeof(Covg))*use_ncols*8 + 1);
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, max_num_kmers,
                                        false, NULL);

  graph_mem = hash_table_mem(kmers_in_hash,bits_per_kmer,NULL);
  bytes_to_str(graph_mem, 1, graph_mem_str);

  if(graph_mem >= args->mem_to_use)
    die("Not enough memory for graph (requires %s)", graph_mem_str);

  // Fringe nodes
  if(dist > 0)
  {
    fringe_mem = args->mem_to_use - graph_mem;
    num_of_fringe_nodes = fringe_mem / (sizeof(dBNode) * 2);

    bytes_to_str(fringe_mem, 1, fringe_mem_str);
    ulong_to_str(num_of_fringe_nodes, num_fringe_nodes_str);

    status("[memory] fringe nodes: %s (%s)\n",
           num_fringe_nodes_str, fringe_mem_str);

    if(num_of_fringe_nodes < 100)
      die("Not enough memory for the graph search (set -m <mem> higher)");
  }

  if(!futil_is_file_writable(out_path))
    die("Cannot write to output file: %s", out_path);

  // Create db_graph
  // multiple colours may be useful later in pulling out multiple colours
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, use_ncols, use_ncols, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity*use_ncols, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity*use_ncols, sizeof(Covg));

  size_t num_words64 = roundup_bits2words64(db_graph.ht.capacity);
  uint64_t *kmer_mask = calloc2(num_words64, sizeof(uint64_t));

  LoadingStats stats;
  loading_stats_init(&stats);

  // Store edge nodes here
  EdgeNodes list0, list1, *listptr0 = &list0, *listptr1 = &list1, *tmplptr;
  edge_nodes_alloc(&list0, &db_graph, kmer_mask, num_of_fringe_nodes);
  edge_nodes_alloc(&list1, &db_graph, kmer_mask, num_of_fringe_nodes);

  //
  // Load graphs
  //
  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  StrBuf intersect_gname;
  strbuf_alloc(&intersect_gname, 1024);

  size_t tmpinto; boolean tmpflatten;
  for(i = 0; i < num_files; i++) {
    tmpinto = files[i].fltr.intocol;
    tmpflatten = files[i].fltr.flatten;

    if(total_cols > db_graph.num_of_cols) {
      // files[i].fltr.intocol = 0;
      file_filter_update_intocol(&files[i].fltr, 0);
      files[i].fltr.flatten = true;
    }

    graph_load(&files[i], gprefs, &stats);
    // files[i].fltr.intocol = tmpinto;
    file_filter_update_intocol(&files[i].fltr, tmpinto);
    files[i].fltr.flatten = tmpflatten;

    for(j = 0; j < files[i].fltr.ncols; j++) {
      col = files[i].fltr.cols[j];
      graph_info_make_intersect(&files[i].hdr.ginfo[col], &intersect_gname);
    }
  }

  hash_table_print_stats(&db_graph.ht);

  char subgraphstr[] = "subgraph:{";
  strbuf_insert(&intersect_gname, 0, subgraphstr, strlen(subgraphstr));
  strbuf_append_char(&intersect_gname, '}');

  // Load sequence and mark in first pass
  read_t r1, r2;
  if(seq_read_alloc(&r1) == NULL || seq_read_alloc(&r2) == NULL)
    die("Out of memory");

  for(i = 0; i < num_seed_files; i++)
    seq_parse_se(seed_files[i], 0, &r1, &r2, store_read_nodes, &list0);

  status("Read in %zu seed kmers\n", list0.stats.kmers_loaded);

  if(dist > 0)
  {
    char tmpstr[100];
    ulong_to_str(dist, tmpstr);
    status("Extending subgraph by %s\n", tmpstr);

    size_t d;
    for(d = 1; d < dist; d++) {
      for(i = 0; i < listptr0->nbuf.len; i++) {
        store_node_neighbours(listptr0->nbuf.data[i].key, listptr1);
      }
      listptr0->nbuf.len = 0;
      SWAP(listptr0, listptr1, tmplptr);
    }
  }

  // DEV: print stats?
  // LoadingStats final_stats = list0.stats;
  // loading_stats_merge(&final_stats, &list1.stats);

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  edge_nodes_dealloc(&list0);
  edge_nodes_dealloc(&list1);

  status("Pruning untouched nodes...\n");

  if(invert) {
    for(i = 0; i < num_words64; i++)
      kmer_mask[i] = ~kmer_mask[i];
  }

  // Remove nodes that were not flagged
  prune_nodes_lacking_flag(&db_graph, kmer_mask);
  free(kmer_mask);

  hash_table_print_stats(&db_graph.ht);

  // Dump nodes that were flagged
  Edges *intersect_edges = NULL;
  boolean kmers_loaded = true;
  boolean colours_loaded = (total_cols <= db_graph.num_of_cols);

  if(!colours_loaded)
  {
    // Need to reload graph colours - therefore construct edge intersection set
    intersect_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));
    for(i = 0; i < db_graph.ht.capacity; i++)
      intersect_edges[i] = db_node_get_edges_union(&db_graph, i);
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
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
