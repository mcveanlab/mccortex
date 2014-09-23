#include "global.h"
#include "subgraph.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "prune_nodes.h"
#include "supernode.h"
#include "loading_stats.h"
#include "util.h"

typedef struct
{
  const dBGraph *const db_graph; // graph we are operating on
  uint8_t *const kmer_mask; // bitset of visited kmers
  const bool grab_supernodes; // grab entire supernodes or just kmers
  dBNodeBuffer nbufs[2], snode_buf;
  LoadingStats stats;
} SubgraphBuilder;

static void subgraph_builder_alloc(SubgraphBuilder *builder,
                                   size_t num_fringe_nodes,
                                   bool grab_supernodes,
                                   uint8_t *kmer_mask,
                                   const dBGraph *graph)
{
  SubgraphBuilder tmp = {.db_graph = graph,
                         .kmer_mask = kmer_mask,
                         .grab_supernodes = grab_supernodes};

  memcpy(builder, &tmp, sizeof(SubgraphBuilder));
  db_node_buf_alloc(&builder->nbufs[0], num_fringe_nodes);
  db_node_buf_alloc(&builder->nbufs[1], num_fringe_nodes);
  db_node_buf_alloc(&builder->snode_buf, 128);
  loading_stats_init(&builder->stats);
}

static void subgraph_builder_dealloc(SubgraphBuilder *builder)
{
  db_node_buf_dealloc(&builder->nbufs[0]);
  db_node_buf_dealloc(&builder->nbufs[1]);
  db_node_buf_dealloc(&builder->snode_buf);
}

// Mark all kmers touched by a read, if they already exist in the graph
static void mark_bkmer(BinaryKmer bkmer, dBNodeBuffer *nbuf,
                       uint8_t *kmer_mask, const dBGraph *db_graph)
{
  dBNode node = db_graph_find(db_graph, bkmer);

  #ifdef CTXVERBOSE
    char tmp[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph->kmer_size, tmp);
    status("got bkmer %s\n", tmp);
  #endif

  // attempt_add return index of item or -1 on failure
  if(node.key != HASH_NOT_FOUND) {
    if(!bitset_get(kmer_mask, node.key) && nbuf->capacity > 0 &&
       db_node_buf_attempt_add(nbuf, node) < 0) {
      die("Please increase <mem> size");
    }
    bitset_set(kmer_mask, node.key);
  }
}

// Mark entire supernodes that are touched by a read
static inline void mark_snode(BinaryKmer bkmer,
                              dBNodeBuffer *nbuf, dBNodeBuffer *snode_buf,
                              uint8_t *kmer_mask, const dBGraph *db_graph)
{
  dBNode node = db_graph_find(db_graph, bkmer);
  size_t i;

  if(node.key != HASH_NOT_FOUND && !bitset_get(kmer_mask, node.key))
  {
    db_node_buf_reset(snode_buf);
    supernode_find(node.key, snode_buf, db_graph);

    for(i = 0; i < snode_buf->len; i++) {
      bitset_set(kmer_mask, snode_buf->data[i].key);
      // attempt_add return index of item or -1 on failure
      if(nbuf->capacity > 0 &&
         db_node_buf_attempt_add(nbuf, snode_buf->data[i]) < 0) {
        die("Please increase <mem> size");
      }
    }
  }
}

static void store_read_nodes(read_t *r1, read_t *r2,
                             uint8_t qoffset1, uint8_t qoffset2, void *ptr)
{
  (void)qoffset1; (void)qoffset2;
  SubgraphBuilder *builder = (SubgraphBuilder*)ptr;
  const dBGraph *db_graph = builder->db_graph;

  if(builder->grab_supernodes)
  {
    READ_TO_BKMERS(r1, db_graph->kmer_size, 0, 0, &builder->stats, mark_snode,
                   &builder->nbufs[0], &builder->snode_buf,
                   builder->kmer_mask, db_graph);
    if(r2 != NULL) {
      READ_TO_BKMERS(r2, db_graph->kmer_size, 0, 0, &builder->stats, mark_snode,
                     &builder->nbufs[0], &builder->snode_buf,
                     builder->kmer_mask, db_graph);
    }
  }
  else
  {
    READ_TO_BKMERS(r1, db_graph->kmer_size, 0, 0, &builder->stats, mark_bkmer,
                   &builder->nbufs[0], builder->kmer_mask, db_graph);
    if(r2 != NULL) {
      READ_TO_BKMERS(r2, db_graph->kmer_size, 0, 0, &builder->stats, mark_bkmer,
                     &builder->nbufs[0], builder->kmer_mask, db_graph);
    }
  }
}

static void store_node_neighbours(const hkey_t hkey, dBNodeBuffer *nbuf,
                                  uint8_t *kmer_mask, const dBGraph *db_graph)
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
  // attempt_add return index of item or -1 on failure
  for(i = 0; i < num_next; i++) {
    if(!bitset_get(kmer_mask, next_nodes[i].key) &&
       db_node_buf_attempt_add(nbuf, next_nodes[i]) < 0) {
      die("Please increase <mem> size");
    }
    bitset_set(kmer_mask, next_nodes[i].key);
  }
}


static void extend(SubgraphBuilder *builder, size_t dist)
{
  const dBGraph *db_graph = builder->db_graph;
  uint8_t *kmer_mask = builder->kmer_mask;
  dBNodeBuffer *nbuf0 = &builder->nbufs[0], *nbuf1 = &builder->nbufs[1];
  size_t d, i;

  if(dist > 0)
  {
    char dist_str[100];
    ulong_to_str(dist, dist_str);
    status("Extending subgraph by %s kmers\n", dist_str);

    for(d = 0; d < dist; d++) {
      db_node_buf_reset(nbuf1);
      for(i = 0; i < nbuf0->len; i++) {
        store_node_neighbours(nbuf0->data[i].key, nbuf1, kmer_mask, db_graph);
      }
      SWAP(nbuf0, nbuf1);
    }
  }
}

static void print_stats(const SubgraphBuilder *builder)
{
  size_t nseed_kmers = builder->stats.num_kmers_loaded;
  size_t nkmers_found = builder->nbufs[0].len;
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

// `nthreads` number of threads to use
// `dist` is how many steps away from seed kmers to take
// `invert`, if true, means only save kmers not touched
// `fringe_mem` is how many bytes can be used to remember the fringe of
// the breadth first search (we use 8 bytes per kmer)
// `kmer_mask` should be a bit array (one bit per kmer) of zero'd memory
void subgraph_from_reads(dBGraph *db_graph, size_t nthreads, size_t dist,
                         bool invert, bool grab_supernodes,
                         size_t fringe_mem, uint8_t *kmer_mask,
                         seq_file_t **files, size_t num_files)
{
  // divide by two since this is the number of nodes per list, two lists
  size_t i, num_of_fringe_nodes = get_num_fringe_nodes(fringe_mem, dist);

  SubgraphBuilder builder;
  subgraph_builder_alloc(&builder, num_of_fringe_nodes, grab_supernodes,
                         kmer_mask, db_graph);

  // Load sequence and mark in first pass
  read_t r1;
  if(seq_read_alloc(&r1) == NULL)
    die("Out of memory");

  for(i = 0; i < num_files; i++)
    seq_parse_se_sf(files[i], 0, &r1, store_read_nodes, &builder);

  print_stats(&builder);

  seq_read_dealloc(&r1);

  extend(&builder, dist);
  subgraph_builder_dealloc(&builder);

  if(invert) {
    status("Inverting selection...");
    size_t num_bytes = roundup_bits2bytes(db_graph->ht.capacity);
    for(i = 0; i < num_bytes; i++)
      kmer_mask[i] = ~kmer_mask[i];
  }

  status("Pruning untouched nodes...");

  prune_nodes_lacking_flag(nthreads, kmer_mask, db_graph);
}

// `nthreads` number of threads to use
// `dist` is how many steps away from seed kmers to take
// `invert`, if true, means only save kmers not touched
// `fringe_mem` is how many bytes can be used to remember the fringe of
// the breadth first search (we use 8 bytes per kmer)
// `kmer_mask` should be a bit array (one bit per kmer) of zero'd memory
void subgraph_from_seq(dBGraph *db_graph, size_t nthreads, size_t dist,
                       bool invert, bool grab_supernodes,
                       size_t fringe_mem, uint8_t *kmer_mask,
                       char **seqs, size_t *seqlens, size_t num_seqs)
{
  // divide by two since this is the number of nodes per list, two lists
  size_t i, num_of_fringe_nodes = get_num_fringe_nodes(fringe_mem, dist);

  SubgraphBuilder builder;
  subgraph_builder_alloc(&builder, num_of_fringe_nodes, grab_supernodes,
                         kmer_mask, db_graph);

  // Load sequence and mark in first pass
  char empty[10] = "";
  read_t r1 = {.name = {.b = empty, .end = 0, .size = 0},
               .seq  = {.b = NULL,  .end = 0, .size = 0},
               .qual = {.b = empty, .end = 0, .size = 0}};

  for(i = 0; i < num_seqs; i++) {
    r1.seq.b = seqs[i];
    r1.seq.end = seqlens[i];
    store_read_nodes(&r1, NULL, 0, 0, &builder);
  }

  print_stats(&builder);

  extend(&builder, dist);
  subgraph_builder_dealloc(&builder);

  if(invert) {
    status("Inverting selection...\n");
    size_t num_bytes = roundup_bits2bytes(db_graph->ht.capacity);
    for(i = 0; i < num_bytes; i++)
      kmer_mask[i] = ~kmer_mask[i];
  }

  status("Pruning untouched nodes...\n");

  prune_nodes_lacking_flag(nthreads, kmer_mask, db_graph);
}
