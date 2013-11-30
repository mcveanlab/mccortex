#include "global.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "supernode.h"


void supernode_reverse(hkey_t *nlist, Orientation *olist, size_t len)
{
  assert(len > 0);
  if(len == 1) { olist[0] = !olist[0]; return; }

  size_t i, j;
  hkey_t tmpnode;
  Orientation tmporient;

  for(i = 0, j=len-1; i <= j; i++, j--) {
    SWAP(nlist[i], nlist[j], tmpnode);
    // swap with reverse has to be done manually
    tmporient = olist[i]; olist[i] = !olist[j]; olist[j] = !tmporient;
  }
}

void supernode_normalise(hkey_t *nlist, Orientation *olist, size_t len)
{
  size_t i;
  assert(len > 0);
  if(nlist[0] > nlist[len-1]) {
    supernode_reverse(nlist, olist, len);
  } else if(nlist[0] == nlist[len-1] && olist[0] == REVERSE) {
    for(i = 0; i < len; i++) olist[i] = !olist[i];
  }
}

// Extend a supernode, nlist[offset] and olist[offset] must already be set
// Walk along nodes starting from node/or, storing the supernode in nlist/olist
// Returns the number of nodes added, adds no more than `limit`
// return -1 if out of space and resize == false
int supernode_extend(const dBGraph *db_graph,
                     hkey_t **nlist, Orientation **olist,
                     size_t offset, size_t *arrlen, boolean resize)
{
  assert(db_graph->num_edge_cols == 1);

  hkey_t node = (*nlist)[offset];
  Orientation orient = (*olist)[offset];
  size_t num_nodes = offset+1, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkey, bkmer = db_graph_oriented_bkmer(db_graph, node, orient);
  const Edges *edges = db_graph->col_edges;

  while(edges_has_precisely_one_edge(edges[node], orient, &nuc))
  {
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);

    bkey = db_node_get_key(bkmer, db_graph->kmer_size);
    node = hash_table_find(&db_graph->ht, bkey);
    orient = db_node_get_orientation(bkey, bkmer);

    assert(node != HASH_NOT_FOUND);

    if(edges_has_precisely_one_edge(edges[node], rev_orient(orient), &nuc))
    {
      if(node == (*nlist)[0] || node == (*nlist)[num_nodes-1]) {
        // don't create a loop A->B->A or a->b->B->A
        break;
      }

      if(num_nodes == *arrlen) {
        if(resize) {
          *arrlen *= 2;
          *nlist = realloc2(*nlist, *arrlen*sizeof(hkey_t));
          *olist = realloc2(*olist, *arrlen*sizeof(Orientation));
        }
        else return -1;
      }

      (*nlist)[num_nodes] = node;
      (*olist)[num_nodes] = orient;
      num_nodes++;
    }
    else break;
  }

  return num_nodes;
}

size_t supernode_find(dBGraph *db_graph, hkey_t node,
                      hkey_t **nlist, Orientation **olist, size_t *arrlen)
{
  int len;
  (*nlist)[0] = node;
  (*olist)[0] = REVERSE;
  len = supernode_extend(db_graph, nlist, olist, 0, arrlen, true);
  supernode_reverse(*nlist, *olist, len);
  len = supernode_extend(db_graph, nlist, olist, len-1, arrlen, true);
  return len;
}

void supernode_print(FILE *out, const dBGraph *db_graph,
                     hkey_t *nodes, Orientation *orients, size_t len)
{
  size_t i, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkmer;
  char tmp[MAX_KMER_SIZE+1];

  bkmer = db_graph_oriented_bkmer(db_graph, nodes[0], orients[0]);
  binary_kmer_to_str(bkmer, kmer_size, tmp);
  fputs(tmp, out);

  for(i = 1; i < len; i++) {
    bkmer = db_node_bkmer(db_graph, nodes[i]);
    nuc = db_node_last_nuc(bkmer, orients[i], kmer_size);
    fputc(binary_nuc_to_char(nuc), out);
  }
}

void supernode_gzprint(gzFile out, const dBGraph *db_graph,
                       hkey_t *nodes, Orientation *orients, size_t len)
{
  size_t i, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkmer;
  char tmp[MAX_KMER_SIZE+1];

  bkmer = db_graph_oriented_bkmer(db_graph, nodes[0], orients[0]);
  binary_kmer_to_str(bkmer, kmer_size, tmp);
  gzputs(out, tmp);

  for(i = 1; i < len; i++) {
    bkmer = db_node_bkmer(db_graph, nodes[i]);
    nuc = db_node_last_nuc(bkmer, orients[i], kmer_size);
    gzputc(out, binary_nuc_to_char(nuc));
  }
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
