#include "global.h"
#include "util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "db_unitig.h"

static bool db_unitig_is_closed_cycle(dBNode n0, BinaryKmer bkey0,
                                      dBNode n1, BinaryKmer bkey1,
                                      const dBGraph *db_graph)
{
  Edges edges0, edges1;
  BinaryKmer shiftkmer;
  Nucleotide nuc;
  const size_t kmer_size = db_graph->kmer_size;

  edges0 = db_node_get_edges_union(db_graph, n0.key);
  if(edges_get_indegree(edges0, n0.orient) != 1) return false;

  edges1 = db_node_get_edges_union(db_graph, n1.key);
  if(edges_get_outdegree(edges1, n1.orient) != 1) return false;

  // Check there is forward edge from last to first
  nuc = bkmer_get_last_nuc(bkey0, n0.orient, kmer_size);
  if(!edges_has_edge(edges1, nuc, n1.orient)) return false;

  shiftkmer = bkmer_shift_add_last_nuc(bkey1, n1.orient, kmer_size, nuc);
  if(binary_kmer_eq(bkey0, shiftkmer)) return true;
  shiftkmer = binary_kmer_reverse_complement(shiftkmer, kmer_size);
  return binary_kmer_eq(bkey0, shiftkmer);
}

// Orient unitig
// Once oriented, unitig has lowest possible kmerkey at the beginning,
// oriented FORWARDs if possible
void db_unitig_normalise(dBNode *nlist, size_t len, const dBGraph *db_graph)
{
  // Sort unitig into forward orientation
  ctx_assert(len > 0);

  if(len == 1) {
    nlist[0].orient = FORWARD;
    return;
  }

  BinaryKmer bkey0 = db_node_get_bkmer(db_graph, nlist[0].key);
  BinaryKmer bkey1 = db_node_get_bkmer(db_graph, nlist[len-1].key);

  // Check if closed cycle
  if(db_unitig_is_closed_cycle(nlist[0], bkey0, nlist[len-1], bkey1, db_graph))
  {
    // find lowest kmer to start from
    BinaryKmer lowest = bkey0, tmp;
    size_t i, lowidx = 0;
    for(i = 1; i < len; i++) {
      tmp = db_node_get_bkmer(db_graph, nlist[i].key);
      if(binary_kmer_lt(tmp, lowest)) {
        lowest = tmp;
        lowidx = i;
      }
    }

    // If already starting from the lowest kmer no change needed
    if(lowidx > 0 || nlist[0].orient != FORWARD)
    {
      // a->b->c->d->e->f->a
      // if c is lowest and FORWARD:  c->d->e->f->a->b (keep orientations)
      // if c is lowest and REVERSE:  c->b->a->f->e->d (reverse orientations)

      if(nlist[lowidx].orient == FORWARD) {
        // Shift left by lowidx, without affecting orientations
        db_nodes_left_shift(nlist, len, lowidx);
      } else {
        db_nodes_reverse_complement(nlist, lowidx+1);
        db_nodes_reverse_complement(nlist+lowidx+1, len-lowidx-1);
      }
    }
  }
  else if(binary_kmer_lt(bkey1, bkey0)) {
    db_nodes_reverse_complement(nlist, len);
  }
}

// Extend a unitig, nlist[offset] must already be set
// Walk along nodes starting from node/or, storing the unitig in nlist
// Returns the number of nodes added, adds no more than `limit`
// return false if out of space and limit > 0
bool db_unitig_extend(dBNodeBuffer *nbuf, size_t limit,
                   const dBGraph *db_graph)
{
  ctx_assert(nbuf->len > 0);

  const size_t kmer_size = db_graph->kmer_size;
  dBNode node0 = nbuf->b[0], node1 = nbuf->b[nbuf->len-1], node = node1;

  BinaryKmer bkmer = db_node_oriented_bkmer(db_graph, node);
  Edges edges = db_node_get_edges_union(db_graph, node.key);
  Nucleotide nuc;

  while(edges_has_precisely_one_edge(edges, node.orient, &nuc))
  {
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    node = db_graph_find(db_graph, bkmer);
    edges = db_node_get_edges_union(db_graph, node.key);

    ctx_assert(node.key != HASH_NOT_FOUND);

    if(edges_has_precisely_one_edge(edges, rev_orient(node.orient), &nuc))
    {
      if(node.key == node0.key || node.key == nbuf->b[nbuf->len-1].key) {
        // don't create a loop A->B->A or a->b->B->A
        break;
      }

      if(limit && nbuf->len >= limit) return false;

      db_node_buf_add(nbuf, node);
    }
    else break;
  }

  return true;
}

void db_unitig_fetch(hkey_t hkey, dBNodeBuffer *nbuf, const dBGraph *db_graph)
{
  dBNode first = {.key = hkey, .orient = REVERSE};
  size_t offset = nbuf->len;
  db_node_buf_add(nbuf, first);
  db_unitig_extend(nbuf, 0, db_graph);
  db_nodes_reverse_complement(nbuf->b+offset, nbuf->len-offset);
  db_unitig_extend(nbuf, 0, db_graph);
}

// Count number of read starts using coverage data
size_t db_unitig_read_starts(const Covg *covgs, size_t len)
{
  if(len == 0) return 0;
  if(len == 1) return covgs[0];

  size_t i, read_starts = covgs[0];

  for(i = 1; i+1 < len; i++)
  {
    if(covgs[i] > covgs[i-1] && covgs[i-1] != covgs[i+1])
      read_starts += covgs[i] - covgs[i-1];
  }

  if(covgs[len-1] > covgs[len-2])
    read_starts += covgs[len-1] - covgs[len-2];

  return read_starts;
}

size_t db_unitig_covg_mean(const Covg *covgs, size_t len)
{
  ctx_assert(len > 0);
  size_t i, sum = 0;
  for(i = 0; i < len; i++) sum += covgs[i];
  return (sum+len/2) / len; // round to nearest integer
}

//
// Iterate over unitigs in the graph with multiple threads
//

static inline int unitig_iterate_node(hkey_t hkey, size_t threadid,
                                      dBNodeBuffer *nbuf,
                                      uint8_t *visited,
                                      const dBGraph *db_graph,
                                      void (*func)(dBNodeBuffer nbuf,
                                                   size_t threadid,
                                                   void *arg),
                                      void *arg)
{
  bool got_lock = false;
  size_t i;

  if(!bitset_get_mt(visited, hkey))
  {
    db_node_buf_reset(nbuf);
    db_unitig_fetch(hkey, nbuf, db_graph);

    // Mark key node (lowest hkey_t value) as visited
    hkey_t node0 = nbuf->b[0].key;
    for(i = 1; i < nbuf->len; i++) node0 = MIN2(node0, nbuf->b[i].key);

    bitlock_try_acquire(visited, node0, &got_lock);

    // Check if someone else marked it first
    if(got_lock)
    {
      // Mark remaining nodes as visited
      // so we don't use them to seed new unitigs
      for(i = 0; i < nbuf->len; i++)
        (void)bitset_set_mt(visited, nbuf->b[i].key);

      func(*nbuf, threadid, arg);
    }
  }

  return 0; // => keep iterating
}

typedef struct {
  const size_t nthreads;
  uint8_t *const visited;
  const dBGraph *db_graph;
  void (*func)(dBNodeBuffer _nbuf, size_t threadid, void *_arg);
  void *arg;
} UnitigIterating;

static void db_unitigs_iterate_thread(void *arg, size_t threadid)
{
  UnitigIterating iter = *(UnitigIterating*)arg;

  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 2048);

  HASH_ITERATE_PART(&iter.db_graph->ht, threadid, iter.nthreads,
                    unitig_iterate_node,
                    threadid, &nbuf, iter.visited, iter.db_graph,
                    iter.func, iter.arg);

  db_node_buf_dealloc(&nbuf);
}

/**
 * @param visited must be initialised to zero, will be dirty upon return
 **/
void db_unitigs_iterate(size_t nthreads, uint8_t *visited,
                        const dBGraph *db_graph,
                        void (*func)(dBNodeBuffer nbuf, size_t threadid, void *arg),
                        void *arg)
{
  UnitigIterating iter = {.nthreads = nthreads,
                          .visited = visited,
                          .db_graph = db_graph,
                          .func = func,
                          .arg = arg};

  util_multi_thread(&iter, nthreads, db_unitigs_iterate_thread);
}
