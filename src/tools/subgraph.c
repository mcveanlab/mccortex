#include "global.h"
#include "subgraph.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "prune_nodes.h"
#include "loading_stats.h"
#include "util.h"

typedef struct
{
  const dBGraph *const db_graph;
  uint64_t *const kmer_mask;
  dBNodeBuffer nbufs[2];
  size_t buf_idx;
  LoadingStats stats;
} EdgeNodes;

static void edge_nodes_alloc(EdgeNodes *enodes, const dBGraph *graph,
                             uint64_t *kmer_mask, size_t capacity)
{
  EdgeNodes tmp = {.db_graph = graph, .kmer_mask = kmer_mask, .buf_idx = 0};
  memcpy(enodes, &tmp, sizeof(EdgeNodes));
  db_node_buf_alloc(&enodes->nbufs[0], capacity);
  db_node_buf_alloc(&enodes->nbufs[1], capacity);
  loading_stats_init(&enodes->stats);
}

static void edge_nodes_dealloc(EdgeNodes *enodes)
{
  db_node_buf_dealloc(&enodes->nbufs[0]);
  db_node_buf_dealloc(&enodes->nbufs[1]);
}

static void mark_bkmer(BinaryKmer bkmer, dBNodeBuffer *nbuf,
                       uint64_t *kmer_mask, const dBGraph *db_graph)
{
  dBNode node = db_graph_find(db_graph, bkmer);

  #ifdef CTXVERBOSE
    char tmp[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph->kmer_size, tmp);
    status("got bkmer %s\n", tmp);
  #endif

  if(node.key != HASH_NOT_FOUND) {
    if(!bitset_get(kmer_mask, node.key) && nbuf->capacity > 0 &&
       !db_node_buf_attempt_add(nbuf, node)) {
      die("Please increase <mem> size");
    }
    bitset_set(kmer_mask, node.key);
  }
}

static void store_read_nodes(read_t *r1, read_t *r2,
                             uint8_t qoffset1, uint8_t qoffset2, void *ptr)
{
  (void)qoffset1; (void)qoffset2;
  EdgeNodes *enodes = (EdgeNodes*)ptr;
  const dBGraph *db_graph = enodes->db_graph;
  dBNodeBuffer *nbuf = &enodes->nbufs[0];

  READ_TO_BKMERS(r1, db_graph->kmer_size, 0, 0, &enodes->stats, mark_bkmer,
                 nbuf, enodes->kmer_mask, db_graph);
  if(r2 != NULL) {
    READ_TO_BKMERS(r2, db_graph->kmer_size, 0, 0, &enodes->stats, mark_bkmer,
                   nbuf, enodes->kmer_mask, db_graph);
  }
}

static void store_node_neighbours(const hkey_t hkey, dBNodeBuffer *nbuf,
                                  uint64_t *kmer_mask, const dBGraph *db_graph)
{
  // Get neighbours
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
    if(!bitset_get(kmer_mask, next_nodes[i].key) &&
       !db_node_buf_attempt_add(nbuf, next_nodes[i])) {
      die("Please increase <mem> size");
    }
    bitset_set(kmer_mask, next_nodes[i].key);
  }
}


static void extend(EdgeNodes *enodes, size_t dist)
{
  const dBGraph *db_graph = enodes->db_graph;
  uint64_t *kmer_mask = enodes->kmer_mask;
  dBNodeBuffer *nbuf0 = &enodes->nbufs[0], *nbuf1 = &enodes->nbufs[1], *tmplist;
  size_t d, i;

  if(dist > 0)
  {
    char tmpstr[100];
    ulong_to_str(dist, tmpstr);
    status("Extending subgraph by %s\n", tmpstr);

    for(d = 0; d < dist; d++) {
      for(i = 0; i < nbuf0->len; i++) {
        store_node_neighbours(nbuf0->data[i].key, nbuf1, kmer_mask, db_graph);
      }
      nbuf0->len = 0;
      SWAP(nbuf0, nbuf1, tmplist);
    }
  }
}

static void print_stats(const EdgeNodes *enodes)
{
  size_t nseed_kmers = enodes->stats.num_kmers_loaded;
  size_t nkmers_found = enodes->nbufs[0].len;
  char nseed_kmers_str[100], nkmers_found_str[100];
  ulong_to_str(nseed_kmers, nseed_kmers_str);
  ulong_to_str(nkmers_found, nkmers_found_str);
  status("Found %s / %s (%.2f%%) seed kmers",
         nkmers_found_str, nseed_kmers_str, (100.0*nkmers_found)/nseed_kmers);
}

static size_t get_num_fringe_nodes(size_t fringe_mem, size_t dist)
{
  return dist == 0 ? 0 : fringe_mem / (sizeof(dBNode) * 2);
}

// `dist` is how many steps away from seed kmers to take
// `invert`, if true, means only save kmers not touched
// `fringe_mem` is how many bytes can be used to remember the fringe of
// the breadth first search (we use 8 bytes per kmer)
// `kmer_mask` should be a bit array (one bit per kmer) of zero'd memory
void subgraph_from_reads(dBGraph *db_graph, size_t dist, bool invert,
                         size_t fringe_mem, uint64_t *kmer_mask,
                         seq_file_t **files, size_t num_files)
{
  // divide by two since this is the number of nodes per list, two lists
  size_t i, num_of_fringe_nodes = get_num_fringe_nodes(fringe_mem, dist);

  EdgeNodes enodes;
  edge_nodes_alloc(&enodes, db_graph, kmer_mask, num_of_fringe_nodes);

  // Load sequence and mark in first pass
  read_t r1;
  if(seq_read_alloc(&r1) == NULL)
    die("Out of memory");

  for(i = 0; i < num_files; i++)
    seq_parse_se_sf(files[i], 0, &r1, store_read_nodes, &enodes);

  print_stats(&enodes);

  seq_read_dealloc(&r1);

  extend(&enodes, dist);
  edge_nodes_dealloc(&enodes);

  if(invert) {
    status("Inverting selection...\n");
    size_t num_words64 = roundup_bits2words64(db_graph->ht.capacity);
    for(i = 0; i < num_words64; i++)
      kmer_mask[i] = ~kmer_mask[i];
  }

  status("Pruning untouched nodes...\n");

  prune_nodes_lacking_flag(db_graph, kmer_mask);
}

// `dist` is how many steps away from seed kmers to take
// `invert`, if true, means only save kmers not touched
// `fringe_mem` is how many bytes can be used to remember the fringe of
// the breadth first search (we use 8 bytes per kmer)
// `kmer_mask` should be a bit array (one bit per kmer) of zero'd memory
void subgraph_from_seq(dBGraph *db_graph, size_t dist, bool invert,
                       size_t fringe_mem, uint64_t *kmer_mask,
                       char **seqs, size_t num_seqs)
{
  // divide by two since this is the number of nodes per list, two lists
  size_t i, num_of_fringe_nodes = get_num_fringe_nodes(fringe_mem, dist);

  EdgeNodes enodes;
  edge_nodes_alloc(&enodes, db_graph, kmer_mask, num_of_fringe_nodes);

  // Load sequence and mark in first pass
  char empty[10] = "";
  read_t r1 = {.name = {.b = empty, .end = 0, .size = 0},
               .seq  = {.b = NULL,  .end = 0, .size = 0},
               .qual = {.b = empty, .end = 0, .size = 0}};

  for(i = 0; i < num_seqs; i++) {
    r1.seq.b = seqs[i];
    r1.seq.end = strlen(seqs[i]);
    store_read_nodes(&r1, NULL, 0, 0, &enodes);
  }

  print_stats(&enodes);

  extend(&enodes, dist);
  edge_nodes_dealloc(&enodes);

  if(invert) {
    status("Inverting selection...\n");
    size_t num_words64 = roundup_bits2words64(db_graph->ht.capacity);
    for(i = 0; i < num_words64; i++)
      kmer_mask[i] = ~kmer_mask[i];
  }

  status("Pruning untouched nodes...\n");

  prune_nodes_lacking_flag(db_graph, kmer_mask);
}
