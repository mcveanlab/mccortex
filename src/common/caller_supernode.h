#ifndef CALLER_SUPERNODE_H_
#define CALLER_SUPERNODE_H_

#include "graph_typedef.h"
#include "hash_table.h"


#define MAX_FLANK_KMERS 1000
#define MAX_ALLELE_KMERS 1000

#define NODE_BUFSIZE(colours) ((MAX_ALLELE_KMERS) * (colours) * 4)

//
// Data types and structs for calling variants with the bubble caller
//
typedef Orientation SuperOrientation;

typedef struct CallerSupernode CallerSupernode;
typedef struct SupernodePath SupernodePath;
typedef struct SupernodePathPos SupernodePathPos;

struct CallerSupernode
{
  hkey_t *nodes;
  Orientation *orients;

  size_t num_of_nodes;

  int num_prev, num_next;
  hkey_t prev_nodes[4], next_nodes[4];
  Orientation prev_orients[4], next_orients[4];

  SupernodePathPos *first_pathpos;
};

struct SupernodePath
{
  int colour;
  CallerSupernode *supernodes[MAX_ALLELE_KMERS];
  SuperOrientation superorients[MAX_ALLELE_KMERS];
  size_t length;
};

struct SupernodePathPos
{
  SupernodePath *path;
  size_t pos;
  SupernodePathPos *next;
};

#define supernode_get_orientation(snode,node,or) \
  ((node) == (snode)->nodes[0] && (or) == (snode)->orients[0] ? forward : reverse)

void reverse_node_list(hkey_t *nlist, Orientation *olist, size_t len);

// Walk along nodes starting from node/or, storing the supernode in nlist/olist
// Returns the number of nodes added, adds no more than `limit`
size_t supernode_traverse(hkey_t node, Orientation or,
                          hkey_t *nlist, Orientation *olist,
                          size_t limit, const dBGraph *db_graph,
                          boolean *out_of_space);

// Create a supernode strating at node/or.  Store in snode. 
// Ensure snode->nodes and snode->orients point to valid memory before passing
// Returns 0 on failure, otherwise snode->num_of_nodes
size_t caller_supernode_create(hkey_t node, Orientation or,
                               CallerSupernode *snode, size_t limit,
                               const dBGraph *db_graph);

#define supernode_pathpos_equal(a,b) (cmp_snpath_pos(a,b) == 0)

uint64_t supernode_pathpos_hash(SupernodePathPos *spp);

// Two SupernodePathPos objects are equal if they describe the same path through
// the graph.  Sort by length (shortest first), then orientation and nodes
int cmp_snpath_pos(const void *p1, const void *p2);


#endif
