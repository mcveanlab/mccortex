#include "global.h"

#include "seq_file.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "graph_format.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" extend [options] <in.ctx> <in.fa> <dist> <out.fa>\n"
"  options:\n"
"   -m <mem>   Memory to use\n"
"   -n <kmers> Kmers in the hash\n";

typedef struct
{
  dBGraph *db_graph;
  dBNodeBuffer *readbuffw, *readbufrv;
  FILE *out;
  char *buf;
  size_t buflen;
  size_t max_walk;
  uint64_t *visited;
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
  BinaryKmer bkey = db_node_get_key(bkmer, db_graph->kmer_size);
  hkey_t node = hash_table_find(&db_graph->ht, bkey);
  if(node == HASH_NOT_FOUND) node_not_in_graph(bkmer, db_graph);
  Orientation orient = db_node_get_orientation(bkmer, bkey);
  nodebuf->data[nodebuf->len].key = node;
  nodebuf->data[nodebuf->len].orient = orient;
  nodebuf->len++;
}

static void walk_graph(hkey_t node, Orientation orient,
                       dBNodeBuffer *buf, size_t max,
                       dBGraph *db_graph, uint64_t *visited)
{
  size_t i, origlen = buf->len;
  Edges edges;
  Nucleotide nuc;
  BinaryKmer bkmer;

  for(i = 0; i < max; i++)
  {
    edges = db_node_edges_union(db_graph, node);
    if(!edges_has_precisely_one_edge(edges, orient, &nuc)) break;
    bkmer = db_node_bkmer(db_graph, node);
    db_graph_next_node(db_graph, bkmer, nuc, orient, &node, &orient);
    if(db_node_has_traversed(visited, node, orient)) break;
    db_node_set_traversed(visited, node, orient);

    buf->data[buf->len].key = node;
    buf->data[buf->len].orient = orient;
    buf->len++;
  }

  for(i = origlen; i < buf->len; i++)
    db_node_fast_clear_traversed(visited, buf->data[i].key);
}

static void extend_read(read_t *r, ExtendContig *contig, SeqLoadingStats *stats)
{
  dBGraph *db_graph = contig->db_graph;
  dBNodeBuffer *readbuffw = contig->readbuffw, *readbufrv = contig->readbufrv;

  if(r->seq.end < db_graph->kmer_size) return;

  db_node_buf_ensure_capacity(readbuffw, r->seq.end + contig->max_walk * 2);
  readbuffw->len = readbufrv->len = 0;
  READ_TO_BKMERS(r, db_graph->kmer_size, 0, 0, stats, add_kmer, db_graph, readbuffw);

  if(readbuffw->len == 0) die("Invalid entry: %s", r->name.b);

  size_t i, j;

  // extend forwards and backwards
  walk_graph(readbuffw->data[readbuffw->len-1].key,
             readbuffw->data[readbuffw->len-1].orient,
             readbuffw, contig->max_walk, db_graph, contig->visited);

  walk_graph(readbuffw->data[0].key,
             opposite_orientation(readbuffw->data[0].orient),
             readbufrv, contig->max_walk, db_graph, contig->visited);

  // Shift forward list up
  memmove(readbuffw->data+readbufrv->len, readbuffw->data,
          readbuffw->len * sizeof(dBNode));

  // Reverse orientation for backwards kmers
  for(i = 0, j = readbufrv->len-1; i < readbufrv->len; i++, j--) {
    readbuffw->data[i].key = readbufrv->data[j].key;
    readbuffw->data[i].orient = opposite_orientation(readbufrv->data[j].orient);
  }

  // to string and print
  size_t len = readbuffw->len + readbufrv->len;
  if(contig->buflen < len+1) {
    contig->buflen = ROUNDUP2POW(len+1);
    contig->buf = realloc2(contig->buf, contig->buflen);
  }

  db_nodes_to_str(readbuffw->data, len, db_graph, contig->buf);

  fprintf(contig->out, ">%s\n%s\n", r->name.b, contig->buf);
}

static void extend_reads(read_t *r1, read_t *r2,
                         int qoffset1, int qoffset2,
                         const SeqLoadingPrefs *prefs, SeqLoadingStats *stats,
                         void *ptr)
{
  (void)qoffset1; (void)qoffset2; (void)prefs;
  ExtendContig *contig = (ExtendContig*)ptr;
  extend_read(r1, contig, stats);
  if(r2 != NULL) extend_read(r2, contig, stats);
}

int ctx_extend(CmdArgs *args)
{
  cmd_accept_options(args, "mn", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc != 4) print_usage(usage, NULL);

  char *input_ctx_path, *input_fa_path, *out_fa_path;
  size_t dist;
  seq_file_t *seq_fa_file;

  input_ctx_path = argv[2];
  input_fa_path = argv[3];

  // Probe graph file
  GraphFileReader file = INIT_GRAPH_READER;
  int ret = graph_file_open(&file, input_ctx_path, false);

  if(ret == 0)
    print_usage(usage, "Cannot read input graph file: %s", input_ctx_path);
  else if(ret < 0)
    print_usage(usage, "Input graph file isn't valid: %s", input_ctx_path);

  if((seq_fa_file = seq_open(input_fa_path)) == NULL)
    print_usage(usage, "Cannot read input FASTA/FASTQ/SAM/BAM file: %s", input_fa_path);

  if(!parse_entire_size(argv[4], &dist))
    print_usage(usage, "Invalid dist argument: %s", argv[4]);

  out_fa_path = argv[5];

  if(!test_file_writable(out_fa_path))
    print_usage(usage, "Cannot write output file: %s", out_fa_path);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;

  bits_per_kmer = (sizeof(Edges) + 2*sizeof(uint64_t)) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        file.hdr.num_of_kmers, true);

  status("Max walk: %zu\n", dist);

  // Set up dBGraph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, file.hdr.kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));

  size_t visited_words = 2*round_bits_to_words64(db_graph.ht.capacity);
  uint64_t *visited = calloc2(visited_words, sizeof(Edges));

    // Store edge nodes here
  dBNodeBuffer readbuffw, readbufrv;
  db_node_buf_alloc(&readbuffw, 2048);
  db_node_buf_alloc(&readbufrv, 1024);

  size_t buflen = 1024;
  char *buf = malloc2(buflen * sizeof(char));

  FILE *out = fopen(out_fa_path, "w");
  if(out == NULL) die("Cannot open output file: %s", out_fa_path);

  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           // binaries
                           .boolean_covgs = false,
                           .must_exist_in_graph = false,
                           .empty_colours = false,
                           // sequence
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false};

  file.fltr.flatten = true;
  file.fltr.intocol = 0;

  // Load binary
  graph_load(&file, &prefs, stats);


  ExtendContig contig = {.db_graph = &db_graph,
                         .readbuffw = &readbuffw, .readbufrv = &readbufrv,
                         .out = out, .buf = buf, .buflen = buflen,
                         .max_walk = dist, .visited = visited};

  // Parse sequence
  read_t r1, r2;
  if(seq_read_alloc(&r1) == NULL || seq_read_alloc(&r2) == NULL)
    die("Out of memory");

  seq_parse_se_sf(seq_fa_file, &r1, &r2, &prefs, stats, extend_reads, &contig);
  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);

  fclose(out);

  free(buf);
  db_node_buf_dealloc(&readbuffw);
  db_node_buf_dealloc(&readbufrv);
  seq_loading_stats_free(stats);
  free(visited);
  free(db_graph.col_edges);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&file);

  return EXIT_SUCCESS;
}
