#include "global.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "binary_kmer.h"
#include "supernode.h"

static size_t extend_supernode(hkey_t init_node, Orientation init_or,
                               hkey_t *end_node, Orientation *end_or,
                               const dBGraph *db_graph, boolean *cycle)
{
  const Edges *edges = db_graph->col_edges;
  Nucleotide nuc;
  hkey_t node = init_node, next_node;
  Orientation or = init_or, next_or;
  size_t len = 0;

  while(edges_has_precisely_one_edge(edges[node], or, &nuc))
  {
    db_graph_next_node_orient(db_graph, db_node_bkmer(db_graph, node), nuc, or,
                              &next_node, &next_or);

    if(!edges_has_precisely_one_edge(edges[next_node], rev_orient(next_or), &nuc))
      break;

    // Check if hit a loop
    if(next_node == init_node) {
      *cycle = true;
      break;
    }

    len++;
    node = next_node;
    or = next_or;
  }

  *end_node = node;
  *end_or = or;

  return len;
}

// In the case of a cycle start node is lowest bkmer
static size_t supernode_cycle(hkey_t init_node, const dBGraph *db_graph,
                              hkey_t *start_node, Orientation *start_orient,
                              hkey_t *end_node, Orientation *end_orient)
{
  const Edges *edges = db_graph->col_edges;

  // Walk forward until we hit init_node
  // store lowest kmer
  hkey_t lowest_node;
  BinaryKmer lowest_bkmer;

  lowest_node = init_node;
  lowest_bkmer = db_node_bkmer(db_graph, init_node);

  hkey_t node = init_node;
  Orientation or = FORWARD;
  Nucleotide nuc;
  size_t len = 0;

  while(edges_has_precisely_one_edge(edges[node], or, &nuc))
  {
    BinaryKmer bkmer = db_node_bkmer(db_graph, node);

    db_graph_next_node_orient(db_graph, bkmer, nuc, or,
                              &node, &or);

    if(node == init_node) break;
    else len++;

    if(binary_kmers_cmp(lowest_bkmer, bkmer) > 0)
    {
      lowest_node = node;
      lowest_bkmer = bkmer;
    }
  }

  // left_node = lowest_node, left_or = FORWARD
  // walk lowest_node in reverse to get right hand node
  *start_node = lowest_node;
  *start_orient = FORWARD;

  edges_has_precisely_one_edge(edges[node], REVERSE, &nuc);
  db_graph_next_node_orient(db_graph, db_node_bkmer(db_graph, lowest_node),
                            nuc, REVERSE, end_node, end_orient);
  *end_orient = rev_orient(*end_orient);

  return len;
}

// Load a supernode from a given node
void supernode_load(hkey_t init_node, const dBGraph *db_graph, Supernode *supernode)
{
  hkey_t left_node, right_node;
  Orientation left_or, right_or;
  size_t len = 1;
  boolean loop = false;

  // Walk forwards
  len += extend_supernode(init_node, FORWARD, &right_node, &right_or,
                          db_graph, &loop);

  if(loop)
  {
    // detected a cycle
    message("Loop detected\n");
    supernode_cycle(init_node, db_graph,
                    &supernode->start_node, &supernode->start_orient,
                    &supernode->end_node, &supernode->end_orient);
  }
  else
  {
    // Walk backwards
    len += extend_supernode(init_node, REVERSE, &left_node, &left_or,
                            db_graph, &loop);

    supernode->start_node = left_node;
    supernode->start_orient = left_or;
    supernode->end_node = right_node;
    supernode->end_orient = right_or;
  }

  supernode->len = len;
}
