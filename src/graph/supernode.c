#include "global.h"
#include "util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "supernode.h"

static bool supernode_is_closed_cycle(const dBNode *nlist, size_t len,
                                         BinaryKmer bkmer0, BinaryKmer bkmer1,
                                         const dBGraph *db_graph)
{
  Edges edges0, edges1;
  BinaryKmer shiftkmer;
  Nucleotide nuc;
  const size_t kmer_size = db_graph->kmer_size;

  edges0 = db_node_get_edges_union(db_graph, nlist[0].key);
  if(edges_get_indegree(edges0, nlist[0].orient) != 1) return false;

  edges1 = db_node_get_edges_union(db_graph, nlist[len-1].key);
  if(edges_get_outdegree(edges1, nlist[len-1].orient) != 1) return false;

  nuc = bkmer_get_last_nuc(bkmer0, nlist[0].orient, kmer_size);
  shiftkmer = bkmer_shift_add_last_nuc(bkmer1, nlist[len-1].orient, kmer_size, nuc);

  if(binary_kmers_are_equal(bkmer0, shiftkmer)) return true;

  shiftkmer = binary_kmer_reverse_complement(shiftkmer, kmer_size);
  return binary_kmers_are_equal(bkmer0, shiftkmer);
}

// Orient supernode
// Once oriented, supernode has lowest possible kmerkey at the beginning,
// oriented FORWARDs if possible
void supernode_normalise(dBNode *nlist, size_t len, const dBGraph *db_graph)
{
  // Sort supernode into forward orientation
  ctx_assert(len > 0);

  if(len == 1) {
    nlist[0].orient = FORWARD;
    return;
  }

  BinaryKmer bkmer0 = db_node_get_bkmer(db_graph, nlist[0].key);
  BinaryKmer bkmer1 = db_node_get_bkmer(db_graph, nlist[len-1].key);

  // Check if closed cycle
  if(supernode_is_closed_cycle(nlist, len, bkmer0, bkmer1, db_graph))
  {
    // find lowest kmer to start from
    BinaryKmer lowest = bkmer0, tmp;
    size_t i, idx = 0;
    for(i = 1; i < len; i++) {
      tmp = db_node_get_bkmer(db_graph, nlist[i].key);
      if(binary_kmer_less_than(tmp, lowest)) {
        lowest = tmp;
        idx = i;
      }
    }

    // If already starting from the lowest kmer no change needed
    if(idx > 0 || nlist[0].orient != FORWARD)
    {
      // a->b->c->d->e->f->a
      // if c is lowest and FORWARD:  c->d->e->f->a->b (keep orientations)
      // if c is lowest and REVERSE:  c->b->a->f->e->d (reverse orientations)

      if(nlist[idx].orient == FORWARD) {
        // Shift left by idx, without affecting orientations
        db_nodes_left_shift(nlist, len, idx);
      } else {
        db_nodes_reverse_complement(nlist, idx+1);
        db_nodes_reverse_complement(nlist+idx+1, len-idx-1);
      }
    }
  }
  else if(binary_kmer_less_than(bkmer1,bkmer0)) {
    db_nodes_reverse_complement(nlist, len);
  }
}

// Extend a supernode, nlist[offset] must already be set
// Walk along nodes starting from node/or, storing the supernode in nlist
// Returns the number of nodes added, adds no more than `limit`
// return false if out of space and limit > 0
bool supernode_extend(dBNodeBuffer *nbuf, size_t limit,
                         const dBGraph *db_graph)
{
  ctx_assert(nbuf->len > 0);

  const size_t kmer_size = db_graph->kmer_size;
  dBNode node0 = nbuf->data[0], node1 = nbuf->data[nbuf->len-1], node = node1;

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
      if(node.key == node0.key || node.key == nbuf->data[nbuf->len-1].key) {
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

void supernode_find(hkey_t hkey, dBNodeBuffer *nbuf, const dBGraph *db_graph)
{
  dBNode first = {.key = hkey, .orient = REVERSE};
  size_t offset = nbuf->len;
  db_node_buf_add(nbuf, first);
  supernode_extend(nbuf, 0, db_graph);
  db_nodes_reverse_complement(nbuf->data+offset, nbuf->len-offset);
  supernode_extend(nbuf, 0, db_graph);
}

// Count number of read starts using coverage data
size_t supernode_read_starts(const Covg *covgs, size_t len)
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

size_t supernode_covg_mean(const Covg *covgs, size_t len)
{
  ctx_assert(len > 0);
  size_t i, sum = 0;
  for(i = 0; i < len; i++) sum += covgs[i];
  return (sum+len/2) / len; // round to nearest integer
}

//
// Iterate over supernodes in the graph with multiple threads
//

static inline int supernode_iterate_node(hkey_t hkey, size_t threadid,
                                         dBNodeBuffer *nbuf,
                                         uint8_t *visited,
                                         const dBGraph *db_graph,
                                         void (*func)(dBNodeBuffer _nbuf,
                                                      size_t threadid,
                                                      void *_arg),
                                         void *arg)
{
  bool got_lock = false;
  size_t i;

  if(!bitset_get_mt(visited, hkey))
  {
    db_node_buf_reset(nbuf);
    supernode_find(hkey, nbuf, db_graph);

    // Mark first node (lowest hkey_t value) as visited
    hkey_t node0 = MIN2(nbuf->data[0].key, nbuf->data[nbuf->len-1].key);
    bitlock_try_acquire(visited, node0, &got_lock);

    // Check if someone else marked it first
    if(got_lock)
    {
      // Mark remaining nodes as visited
      for(i = 0; i < nbuf->len; i++)
        (void)bitset_set_mt(visited, nbuf->data[i].key);

      func(*nbuf, threadid, arg);
    }
  }

  return 0; // => keep iterating
}

typedef struct {
  const size_t threadid, nthreads;
  uint8_t *const visited;
  const dBGraph *db_graph;
  void (*func)(dBNodeBuffer _nbuf, size_t threadid, void *_arg);
  void *arg;
} SupernodeIterator;

static void supernodes_iterate_thread(void *arg)
{
  SupernodeIterator cl = *(SupernodeIterator*)arg;

  dBNodeBuffer nbuf;
  db_node_buf_alloc(&nbuf, 2048);

  HASH_ITERATE_PART(&cl.db_graph->ht, cl.threadid, cl.nthreads,
                    supernode_iterate_node,
                    cl.threadid, &nbuf, cl.visited, cl.db_graph,
                    cl.func, cl.arg);

  db_node_buf_dealloc(&nbuf);
}

void supernodes_iterate(size_t nthreads, uint8_t *visited,
                        const dBGraph *db_graph,
                        void (*func)(dBNodeBuffer _nbuf,
                                     size_t threadid,
                                     void *_arg),
                        void *arg)
{
  size_t i;
  SupernodeIterator *workers = ctx_calloc(nthreads, sizeof(SupernodeIterator));

  for(i = 0; i < nthreads; i++) {
    SupernodeIterator tmp = {.threadid = i, .nthreads = nthreads,
                             .visited = visited, .db_graph = db_graph,
                             .func = func, .arg = arg};
    memcpy(&workers[i], &tmp, sizeof(SupernodeIterator));
  }

  util_run_threads(workers, nthreads, sizeof(SupernodeIterator),
                   nthreads, supernodes_iterate_thread);

  ctx_free(workers);
}
