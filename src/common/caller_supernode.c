#include "global.h"
#include "caller_supernode.h"
#include "db_graph.h"
#include "db_node.h"

#ifdef DEBUG
#define DEBUG_CALLER 1
#endif

void reverse_node_list(hkey_t *nlist, Orientation *olist, size_t len)
{
  if(len == 0) return;
  else if(len == 1) olist[0] = opposite_orientation(olist[0]);
  else
  {
    size_t i, j;
    hkey_t tmp_n;
    Orientation tmp_o;
    for(i = 0, j = len-1; i <= j; i++, j--)
    {
      SWAP(nlist[i], nlist[j], tmp_n);
      tmp_o = olist[i];
      olist[i] = opposite_orientation(olist[j]);
      olist[j] = opposite_orientation(tmp_o);
    }
  }
}

// Orient supernode
static void supernode_naturalise(hkey_t *nlist, Orientation *olist, size_t len)
{
  // Sort supernode into forward orientation
  if(len == 1)
    olist[0] = FORWARD;
  else if(nlist[0] > nlist[len-1] ||
          (nlist[0] == nlist[len-1] && olist[0] > olist[len-1]))
  {
    reverse_node_list(nlist, olist, len);
  }
}

// Walk along nodes starting from node/or, storing the supernode in nlist/olist
// Returns the number of nodes added, adds no more than `limit`
size_t supernode_traverse(hkey_t node, Orientation orient,
                          hkey_t *nlist, Orientation *olist,
                          size_t limit, const dBGraph *db_graph,
                          boolean *out_of_space)
{
  size_t num_nodes = 1, kmer_size = db_graph->kmer_size;
  Nucleotide nuc;
  BinaryKmer bkey, bkmer = db_graph_oriented_bkmer(db_graph, node, orient);
  const Edges *edges = db_graph->edges;

  #ifdef DEBUG_CALLER
    char tmpstr1[MAX_KMER_SIZE+1], tmpstr2[MAX_KMER_SIZE+1];
    binary_kmer_to_str(bkmer, db_graph->kmer_size, tmpstr1);
    printf("  fetch %s:%i\n", tmpstr1, (int)orient);
  #endif

  nlist[0] = node;
  olist[0] = orient;

  while(edges_has_precisely_one_edge(edges[node], orient, &nuc))
  {
    binary_kmer_left_shift_add(&bkmer, kmer_size, nuc);

    bkey = db_node_get_key(bkmer, db_graph->kmer_size);
    node = hash_table_find(&db_graph->ht, bkey);
    orient = db_node_get_orientation(bkey, bkmer);

    #ifdef DEBUG_CALLER
      binary_kmer_to_str(bkmer, db_graph->kmer_size, tmpstr1);
      binary_kmer_to_str(bkey, db_graph->kmer_size, tmpstr2);
      printf(">%s %s:%i nuc:%c\n", tmpstr1, tmpstr2, orient, binary_nuc_to_char(nuc));
    #endif

    assert(node != HASH_NOT_FOUND);

    if(edges_has_precisely_one_edge(edges[node], rev_orient(orient), &nuc))
    {
      if(num_nodes == limit) { *out_of_space = 1; break; }
      if(node == nlist[0] && orient == olist[0]) break; // don't create a loop

      nlist[num_nodes] = node;
      olist[num_nodes] = orient;
      num_nodes++;
    }
    else break;
  }

  // char tmp[1000];
  // nodes_to_str(nlist, olist, num_nodes, db_graph->kmer_size, tmp);
  // printf("   supernode: %s\n", tmp);

  return num_nodes;
}

// Create a supernode strating at node/or.  Store in snode.
// Ensure snode->nodes and snode->orients point to valid memory before passing
// Returns 0 on failure, otherwise snode->num_of_nodes
size_t caller_supernode_create(hkey_t node, Orientation or,
                               CallerSupernode *snode, size_t limit,
                               const dBGraph *db_graph)
{
  Nucleotide nuc;
  Edges union_edges;
  hkey_t first_node, last_node;
  Orientation first_or, last_or;
  BinaryKmer bkmer;

  #ifdef DEBUG_CALLER
    char tmpstr[MAX_KMER_SIZE+1];
    bkmer = db_node_bkmer(db_graph, node);
    binary_kmer_to_str(bkmer, db_graph->kmer_size, tmpstr);
    printf(" create %s:%i\n", tmpstr, (int)or);
  #endif

  // extend path
  boolean out_of_space = 0;
  snode->num_of_nodes = supernode_traverse(node, or, snode->nodes, snode->orients,
                                           limit, db_graph, &out_of_space);

  if(out_of_space)
    return 0;

  supernode_naturalise(snode->nodes, snode->orients, snode->num_of_nodes);

  snode->first_pathpos = NULL;
  snode->num_prev = 0;
  snode->num_next = 0;

  first_node = snode->nodes[0];
  first_or = opposite_orientation(snode->orients[0]);
  last_node = snode->nodes[snode->num_of_nodes-1];
  last_or = snode->orients[snode->num_of_nodes-1];

  // Prev nodes
  union_edges = edges_with_orientation(db_graph->edges[first_node], first_or);
  bkmer = db_node_bkmer(db_graph, first_node);

  for(nuc = 0; nuc < 4; nuc++) {
    if(edges_has_edge(union_edges, nuc, FORWARD)) {
      db_graph_next_node(db_graph, bkmer, nuc, first_or,
                         snode->prev_nodes + snode->num_prev,
                         snode->prev_orients + snode->num_prev);
      snode->num_prev++;
    }
  }

  // Next nodes
  union_edges = edges_with_orientation(db_graph->edges[last_node], last_or);
  bkmer = db_node_bkmer(db_graph, last_node);

  for(nuc = 0; nuc < 4; nuc++) {
    if(edges_has_edge(union_edges, nuc, FORWARD)) {
      db_graph_next_node(db_graph, bkmer, nuc, last_or,
                         snode->next_nodes + snode->num_next,
                         snode->next_orients + snode->num_next);
      snode->num_next++;
    }
  }

  #ifdef DEBUG_CALLER
    char tmpstr1[MAX_KMER_SIZE+1], tmpstr2[MAX_KMER_SIZE+1];
    BinaryKmer first_bkmer = db_node_bkmer(db_graph, first_node);
    BinaryKmer last_bkmer = db_node_bkmer(db_graph, last_node);
    binary_kmer_to_str(first_bkmer, db_graph->kmer_size, tmpstr1);
    binary_kmer_to_str(last_bkmer, db_graph->kmer_size, tmpstr2);
    printf("   ( [>%i] first:%s:%i; len:%zu last:%s:%i [%i<] )\n",
           snode->num_prev, tmpstr1, (int)first_or, snode->num_of_nodes,
           tmpstr2, (int)last_or, snode->num_next);
  #endif

  return snode->num_of_nodes;
}

uint64_t supernode_pathpos_hash(SupernodePathPos *spp)
{
  uint32_t hsh, len = spp->pos + 1;
  size_t snode_size = sizeof(CallerSupernode*) * len;
  size_t sorients_size = sizeof(SuperOrientation) * len;

#ifdef CITY_HASH
  // Use Google's CityHash
  hsh = CityHash32((char*)spp->path->supernodes, snode_size);
  hsh ^= CityHash32((char*)spp->path->superorients, sorients_size);
#else
  // Use Bob Jenkin's lookup3
  hsh = hashlittle(spp->path->supernodes, snode_size, 0);
  hsh = hashlittle(spp->path->superorients, sorients_size, hsh);
#endif

  return hsh;
}

// Two SupernodePathPos objects are equal if they describe the same path through
// the graph.  Sort by length (shortest first), then orientation and nodes
int cmp_snpath_pos(const void *p1, const void *p2)
{
  const SupernodePathPos * const pp1 = p1;
  const SupernodePathPos * const pp2 = p2;

  ptrdiff_t cmp = (ptrdiff_t)pp1->pos - pp2->pos;

  if(cmp != 0)
    return cmp;

  int i, len = pp1->pos;
  for(i = 0; i < len; i++)
  {
    cmp = pp1->path->superorients[i] - pp2->path->superorients[i];

    if(cmp != 0)
      return cmp;

    cmp = pp1->path->supernodes[i] - pp2->path->supernodes[i];

    if(cmp != 0)
      return (cmp > 0 ? 1 : -1);
  }

  return 0;
}
