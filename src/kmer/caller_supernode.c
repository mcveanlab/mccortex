#include "global.h"
#include "caller_supernode.h"
#include "db_graph.h"
#include "db_node.h"
#include "supernode.h"

#ifdef CITY_HASH
  #include "city.h"
#else
  #include "lookup3.h"
#endif

#ifdef CTXVERBOSE
  #define DEBUG_CALLER 1
#endif

// Create a supernode starting at node/or.  Store in snode.
// Ensure snode->nodes and snode->orients point to valid memory before passing
// Returns 0 on failure, otherwise snode->num_of_nodes
size_t caller_supernode_create(dBNode node, CallerSupernode *snode,
                               const dBGraph *db_graph)
{
  assert(db_graph->num_edge_cols == 1);

  #ifdef DEBUG_CALLER
    char tmpstr[MAX_KMER_SIZE+1];
    bkmer = db_node_get_bkmer(db_graph, node.key);
    binary_kmer_to_str(bkmer, db_graph->kmer_size, tmpstr);
    printf(" create %s:%i\n", tmpstr, (int)orient);
  #endif

  dBNodeBuffer *nbuf = snode->nbuf;
  snode->nbuf_offset = nbuf->len;
  db_node_buf_add(nbuf, node);
  supernode_extend(nbuf, 0, db_graph);
  snode->num_of_nodes = nbuf->len - snode->nbuf_offset;
  snode->first_pathpos = NULL;

  dBNode *nodes = snode_nodes(snode);
  supernode_normalise(nodes, snode->num_of_nodes);

  Edges union_edges;
  BinaryKmer bkmer;
  dBNode first, last;

  first = db_node_reverse(nodes[0]);
  last = nodes[snode->num_of_nodes-1];

  // prev nodes
  union_edges = db_node_get_edges_union(db_graph, first.key);
  bkmer = db_node_get_bkmer(db_graph, first.key);
  snode->num_prev = db_graph_next_nodes(db_graph, bkmer, first.orient, union_edges,
                                        snode->prev_nodes, snode->prev_bases);

  // next nodes
  union_edges = db_node_get_edges_union(db_graph, last.key);
  bkmer = db_node_get_bkmer(db_graph, last.key);
  snode->num_next = db_graph_next_nodes(db_graph, bkmer, last.orient, union_edges,
                                        snode->next_nodes, snode->next_bases);

  #ifdef DEBUG_CALLER
    char tmpstr1[MAX_KMER_SIZE+1], tmpstr2[MAX_KMER_SIZE+1];
    BinaryKmer first_bkmer = db_node_get_bkmer(db_graph, first.key);
    BinaryKmer last_bkmer = db_node_get_bkmer(db_graph, last.key);
    binary_kmer_to_str(first_bkmer, db_graph->kmer_size, tmpstr1);
    binary_kmer_to_str(last_bkmer, db_graph->kmer_size, tmpstr2);
    printf("   ( [>%i] first:%s:%i; len:%zu last:%s:%i [%i<] )\n",
           snode->num_prev, tmpstr1, (int)first.orient, snode->num_of_nodes,
           tmpstr2, (int)last.orient, snode->num_next);
  #endif

  return snode->num_of_nodes;
}

uint32_t supernode_pathpos_hash(SupernodePathPos *spp)
{
  uint32_t hsh;
  size_t len = spp->pos + 1;
  size_t snode_size = sizeof(CallerSupernode*) * len;
  size_t sorients_size = sizeof(SuperOrientation) * len;

#ifdef CITY_HASH
  // Use Google's CityHash
  hsh = CityHash32((char*)spp->path->supernodes, snode_size);
  hsh ^= CityHash32((char*)spp->path->superorients, sorients_size);
#else
  // Use Bob Jenkin's lookup3
  hsh = lk3_hashlittle(spp->path->supernodes, snode_size, 0);
  hsh = lk3_hashlittle(spp->path->superorients, sorients_size, hsh);
#endif

  return hsh;
}

// Two SupernodePathPos objects are equal if they describe the same path through
// the graph.  Sort by length (shortest first), then orientation and nodes
int cmp_snpath_pos(const void *p1, const void *p2)
{
  const SupernodePathPos * const pp1 = p1;
  const SupernodePathPos * const pp2 = p2;

  ptrdiff_t cmp = (ptrdiff_t)pp1->pos - (ptrdiff_t)pp2->pos;

  if(cmp != 0)
    return (cmp > 0 ? 1 : -1);

  size_t i, len = pp1->pos;
  for(i = 0; i < len; i++)
  {
    cmp = pp1->path->superorients[i] - pp2->path->superorients[i];

    if(cmp != 0)
      return (cmp > 0 ? 1 : -1);

    cmp = pp1->path->supernodes[i] - pp2->path->supernodes[i];

    if(cmp != 0)
      return (cmp > 0 ? 1 : -1);
  }

  return 0;
}
