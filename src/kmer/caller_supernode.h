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
  dBNodeBuffer *nbuf;
  size_t nbuf_offset, num_of_nodes; // Offset and lenth in nbuf

  // DEV: could save space by using dBNode instead of separate hkey_t+orient
  uint8_t num_prev, num_next;
  hkey_t prev_nodes[4], next_nodes[4];
  Orientation prev_orients[4], next_orients[4];

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

#define supernode_get_orientation(snode,node,or) \
        ((node) == snode_nodes(snode)[0].key && \
         (or) == snode_nodes(snode)[0].orient ? FORWARD : REVERSE)

// Create a supernode strating at node/or.  Store in snode.
size_t caller_supernode_create(hkey_t node, Orientation orient,
                               CallerSupernode *snode, const dBGraph *db_graph);

#define supernode_pathpos_equal(a,b) (cmp_snpath_pos(a,b) == 0)

uint32_t supernode_pathpos_hash(SupernodePathPos *spp);

// Two SupernodePathPos objects are equal if they describe the same path through
// the graph.  Sort by length (shortest first), then orientation and nodes
int cmp_snpath_pos(const void *p1, const void *p2);


#endif
