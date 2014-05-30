#include "global.h"

#include "seq_file.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "seq_reader.h"

const char extend_usage[] =
"usage: "CMD" extend [options] <in.ctx> <in.fa> <dist> <out.fa>\n"
"\n"
"  Extend contigs along supernodes using a populationg graph.\n"
"\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"\n";

typedef struct
{
  dBGraph *db_graph;
  dBNodeBuffer *readbuffw, *readbufrv;
  FILE *out;
  char *buf;
  size_t buflen;
  size_t max_walk;
  uint8_t *visited;
  LoadingStats *stats;
} ExtendContig;

static void node_not_in_graph(BinaryKmer bkmer, dBGraph *db_graph)
__attribute__((noreturn));

static void node_not_in_graph(BinaryKmer bkmer, dBGraph *db_graph)
{
  char str[MAX_KMER_SIZE+1];
  binary_kmer_to_str(bkmer, db_graph->kmer_size, str);
  die("Node not in graph: %s", str);
}

static void add_kmer(BinaryKmer bkmer, dBGraph *db_graph, dBNodeBuffer *nodebuf)
{
  BinaryKmer bkey = bkmer_get_key(bkmer, db_graph->kmer_size);
  hkey_t node = hash_table_find(&db_graph->ht, bkey);
  if(node == HASH_NOT_FOUND) node_not_in_graph(bkmer, db_graph);
  Orientation orient = bkmer_get_orientation(bkmer, bkey);
  nodebuf->data[nodebuf->len].key = node;
  nodebuf->data[nodebuf->len].orient = orient;
  nodebuf->len++;
}

static void walk_graph(dBNode node, dBNodeBuffer *buf, size_t max,
                       dBGraph *db_graph, uint8_t *visited)
{
  size_t i, origlen = buf->len;
  Edges edges;
  Nucleotide nuc;
  BinaryKmer bkmer;

  for(i = 0; i < max; i++)
  {
    edges = db_node_get_edges_union(db_graph, node.key);
    if(!edges_has_precisely_one_edge(edges, node.orient, &nuc)) break;
    bkmer = db_node_get_bkmer(db_graph, node.key);
    node = db_graph_next_node(db_graph, bkmer, nuc, node.orient);

    if(db_node_has_traversed(visited, node)) break;
    db_node_set_traversed(visited, node);

    buf->data[buf->len++] = node;
  }

  for(i = origlen; i < buf->len; i++)
    db_node_fast_clear_traversed(visited, buf->data[i].key);
}

static void extend_read(read_t *r, ExtendContig *contig)
{
  dBGraph *db_graph = contig->db_graph;
  dBNodeBuffer *readbuffw = contig->readbuffw, *readbufrv = contig->readbufrv;

  if(r->seq.end < db_graph->kmer_size) return;

  db_node_buf_ensure_capacity(readbuffw, r->seq.end + contig->max_walk * 2);
  db_node_buf_ensure_capacity(readbufrv, r->seq.end + contig->max_walk);

  readbuffw->len = readbufrv->len = 0;

  READ_TO_BKMERS(r, db_graph->kmer_size, 0, 0, contig->stats,
                 add_kmer, db_graph, readbuffw);

  if(readbuffw->len == 0) die("Invalid entry: %s", r->name.b);

  size_t i, j;

  // extend forwards and backwards
  walk_graph(readbuffw->data[readbuffw->len-1],
             readbuffw, contig->max_walk, db_graph, contig->visited);

  walk_graph(db_node_reverse(readbuffw->data[0]),
             readbufrv, contig->max_walk, db_graph, contig->visited);

  // Shift forward list up
  memmove(readbuffw->data+readbufrv->len, readbuffw->data,
          readbuffw->len * sizeof(dBNode));

  // Reverse orientation for backwards kmers
  for(i = 0, j = readbufrv->len-1; i < readbufrv->len; i++, j--) {
    readbuffw->data[i] = db_node_reverse(readbufrv->data[j]);
    // readbuffw->data[i].key = readbufrv->data[j].key;
    // readbuffw->data[i].orient = opposite_orientation(readbufrv->data[j].orient);
  }

  // to string and print
  size_t len = readbuffw->len + readbufrv->len;
  if(contig->buflen < len+1) {
    contig->buflen = roundup2pow(len+1);
    contig->buf = ctx_realloc(contig->buf, contig->buflen);
  }

  db_nodes_to_str(readbuffw->data, len, db_graph, contig->buf);

  fprintf(contig->out, ">%s\n%s\n", r->name.b, contig->buf);
}

static void extend_reads(read_t *r1, read_t *r2,
                         uint8_t qoffset1, uint8_t qoffset2,
                         void *ptr)
{
  (void)qoffset1; (void)qoffset2;
  ExtendContig *contig = (ExtendContig*)ptr;
  extend_read(r1, contig);
  if(r2 != NULL) extend_read(r2, contig);
}

int ctx_extend(CmdArgs *args)
{
  char **argv = args->argv;
  // Already checked that we have exactly 4 arguments

  char *input_ctx_path, *input_fa_path, *out_fa_path;
  size_t dist;
  seq_file_t *seq_fa_file;

  input_ctx_path = argv[2];
  input_fa_path = argv[3];

  // Probe graph file
  GraphFileReader file = INIT_GRAPH_READER;
  graph_file_open(&file, input_ctx_path, true);

  if((seq_fa_file = seq_open(input_fa_path)) == NULL)
    cmd_print_usage("Cannot read input FASTA/FASTQ/SAM/BAM file: %s", input_fa_path);

  if(!parse_entire_size(argv[4], &dist))
    cmd_print_usage("Invalid dist argument: %s", argv[4]);

  out_fa_path = argv[5];

  if(!futil_is_file_writable(out_fa_path))
    cmd_print_usage("Cannot write output file: %s", out_fa_path);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  // two bits for remembering which orientation we have traversed kmers in
  bits_per_kmer = sizeof(Edges) * 8 + 2;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        file.num_of_kmers, file.num_of_kmers,
                                        true, &graph_mem);

  cmd_check_mem_limit(args->mem_to_use, graph_mem);

  status("Max walk: %zu\n", dist);

  // Set up dBGraph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, file.hdr.kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(Edges));

  uint8_t *visited = ctx_calloc(2*roundup_bits2bytes(db_graph.ht.capacity), 1);

    // Store edge nodes here
  dBNodeBuffer readbuffw, readbufrv;
  db_node_buf_alloc(&readbuffw, 2048);
  db_node_buf_alloc(&readbufrv, 1024);

  size_t buflen = 1024;
  char *buf = ctx_malloc(buflen * sizeof(char));

  FILE *out = fopen(out_fa_path, "w");
  if(out == NULL) die("Cannot open output file: %s", out_fa_path);

  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .empty_colours = false};

  file.fltr.flatten = true;
  file_filter_update_intocol(&file.fltr, 0);

  graph_load(&file, gprefs, &stats);

  ExtendContig contig = {.db_graph = &db_graph,
                         .readbuffw = &readbuffw, .readbufrv = &readbufrv,
                         .out = out, .buf = buf, .buflen = buflen,
                         .max_walk = dist, .visited = visited,
                         .stats = &stats};

  // Parse sequence
  read_t r1;
  if(seq_read_alloc(&r1) == NULL)
    die("Out of memory");

  seq_parse_se_sf(seq_fa_file, 0, &r1, extend_reads, &contig);
  seq_read_dealloc(&r1);
  seq_close(seq_fa_file);

  fclose(out);

  ctx_free(buf);
  db_node_buf_dealloc(&readbuffw);
  db_node_buf_dealloc(&readbufrv);
  ctx_free(visited);
  ctx_free(db_graph.col_edges);
  db_graph_dealloc(&db_graph);

  graph_file_close(&file);

  return EXIT_SUCCESS;
}
