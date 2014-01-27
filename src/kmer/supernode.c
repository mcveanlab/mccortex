#include "global.h"
#include "hash_table.h"
#include "db_graph.h"
#include "binary_kmer.h"
#include "supernode.h"


void supernode_reverse(dBNode *nlist, size_t len)
{
  assert(len > 0);
  size_t i, j;
  dBNode tmpnode;

  for(i = 0, j = len-1; i < j; i++, j--) {
    SWAP(nlist[i], nlist[j], tmpnode);
    nlist[i].orient = !nlist[i].orient;
    nlist[j].orient = !nlist[j].orient;
  }

  if(i == j) nlist[i].orient = !nlist[i].orient;
}

// Orient supernode
void supernode_normalise(dBNode *nlist, size_t len)
{
  // Sort supernode into forward orientation
  assert(len > 0);
  if(len == 1)
    nlist[0].orient = FORWARD;
  else if(nlist[0].key > nlist[len-1].key)
    supernode_reverse(nlist, len);
}

// Extend a supernode, nlist[offset] must already be set
// Walk along nodes starting from node/or, storing the supernode in nlist
// Returns the number of nodes added, adds no more than `limit`
// return false if out of space and limit > 0
boolean supernode_extend(dBNodeBuffer *nbuf, size_t limit, const dBGraph *db_graph)
{
  assert(db_graph->num_edge_cols == 1);

  hkey_t hkey = nbuf->data[nbuf->len-1].key;
  Orientation orient = nbuf->data[nbuf->len-1].orient;
  const size_t kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkey, bkmer = db_graph_oriented_bkmer(db_graph, hkey, orient);
  const Edges *edges = db_graph->col_edges;

  while(edges_has_precisely_one_edge(edges[hkey], orient, &nuc))
  {
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    bkey = db_node_get_key(bkmer, db_graph->kmer_size);
    hkey = hash_table_find(&db_graph->ht, bkey);
    orient = db_node_get_orientation(bkey, bkmer);

    assert(hkey != HASH_NOT_FOUND);

    if(edges_has_precisely_one_edge(edges[hkey], rev_orient(orient), &nuc))
    {
      if(hkey == nbuf->data[0].key || hkey == nbuf->data[nbuf->len-1].key) {
        // don't create a loop A->B->A or a->b->B->A
        break;
      }

      if(limit && nbuf->len >= limit) return false;

      dBNode next = {.key = hkey, .orient = orient};
      db_node_buf_add(nbuf, next);
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
  supernode_reverse(nbuf->data+offset, nbuf->len-offset);
  supernode_extend(nbuf, 0, db_graph);
}

uint32_t supernode_read_starts(uint32_t *covgs, uint32_t len)
{
  if(len == 0) return 0;
  if(len == 1) return covgs[0];

  uint32_t i, read_starts = covgs[0];

  for(i = 1; i+1 < len; i++)
  {
    if(covgs[i] > covgs[i-1] && covgs[i-1] != covgs[i+1])
      read_starts += covgs[i] - covgs[i-1];
  }

  if(covgs[len-1] > covgs[len-2])
    read_starts += covgs[len-1] - covgs[len-2];

  return read_starts;
}
