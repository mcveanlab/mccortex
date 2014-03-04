#ifndef CALLER_SUPERNODE_H_
#define CALLER_SUPERNODE_H_

#include "cortex_types.h"
#include "hash_table.h"
#include "db_node.h"

//
// Data types and structs for calling variants with the bubble caller
//
typedef Orientation SuperOrientation;

typedef struct CallerSupernode CallerSupernode;
typedef struct SupernodePath SupernodePath;
typedef struct SupernodePathPos SupernodePathPos;

struct CallerSupernode
{
  dBNodeBuffer *nbuf; // shared node buffer
  size_t nbuf_offset, num_of_nodes; // Offset and lenth in nbuf
  // Edges to/from this supernode
  dBNode prev_nodes[4], next_nodes[4];
  Nucleotide prev_bases[4], next_bases[4];
  uint8_t num_prev, num_next;
  // Linked list of paths that use this supernode
  SupernodePathPos *first_pathpos;
};

struct SupernodePath
{
  int colour;
  CallerSupernode **supernodes;
  SuperOrientation *superorients;
  size_t length;
};

struct SupernodePathPos
{
  SupernodePath *path;
  size_t pos;
  SupernodePathPos *next;
};

#define snode_nodes(sn)   ((sn)->nbuf->data+(sn)->nbuf_offset)

#define supernode_get_orientation(snode,node) \
        (db_nodes_match(node, snode_nodes(snode)[0]) ? FORWARD : REVERSE)

// Create a supernode starting at node/or.  Store in snode.
// Ensure snode->nodes and snode->orients point to valid memory before passing
// Returns 0 on failure, otherwise snode->num_of_nodes
size_t caller_supernode_create(dBNode node, CallerSupernode *snode,
                               const dBGraph *db_graph);

#define supernode_pathpos_equal(a,b) (cmp_snpath_pos(a,b) == 0)

uint32_t supernode_pathpos_hash(SupernodePathPos *spp);

// Two SupernodePathPos objects are equal if they describe the same path through
// the graph.  Sort by length (shortest first), then orientation and nodes
int cmp_snpath_pos(const void *p1, const void *p2);

#endif
