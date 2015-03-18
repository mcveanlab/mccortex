#include "global.h"
#include "pop_bubbles.h"
#include "supernode.h"

typedef struct
{
  uint8_t *const visited, *const remove;
  dBNodeBuffer *alts;
  const size_t min_covg;
  const dBGraph *db_graph;
} PopBubbles;

// Go forward one node and back one node to get all 'parallel nodes'
// Given node, return nodes {a,b}
//      a
//       \
// node -> x
//       /
//      b
static inline uint8_t get_parallel_nodes(const dBGraph *db_graph, dBNode node,
                                         dBNode nodes[16])
{
  dBNode tmp_nodes[4];
  Nucleotide tmp_nucs[4], tmp_nucs2[4];
  uint8_t i, num_next, tmp_n, n = 0;

  Nucleotide lost_nuc = db_node_get_first_nuc(node, db_graph);
  num_next = db_graph_next_nodes_union(db_graph, node, tmp_nodes, tmp_nucs);

  for(i = 0; i < num_next; i++) {
    tmp_n = db_graph_prev_nodes_with_mask(db_graph, tmp_nodes[i],
                                          lost_nuc, -1,
                                          nodes+n, tmp_nucs2);
    n += tmp_n;
  }

  return n;
}

static inline void mark_node_bitarr_mt(const dBNode *nodes, size_t n,
                                       uint8_t *bitset)
{
  size_t i;
  for(i = 0; i < n; i++) bitset_set_mt(bitset, nodes[i].key);
}

/**
 * Remove the lowest mean coverage branch, by marking the remove bit array.
 * @param min_covg keep all branches with mean coverage >= min_covg
 */
static inline void found_bubble(const dBNode *s1, size_t n1,
                                const dBNode *s2, size_t n2,
                                size_t min_covg,
                                uint8_t *visited, uint8_t *remove,
                                const dBGraph *db_graph)
{
  size_t i, sum_covg1 = 0, sum_covg2 = 0;

  for(i = 0; i < n1; i++) sum_covg1 += db_node_sum_covg(db_graph, s1[i].key);
  for(i = 0; i < n2; i++) sum_covg2 += db_node_sum_covg(db_graph, s2[i].key);
  sum_covg1 /= n1;
  sum_covg2 /= n2;

  // Print alleles
  // pthread_mutex_lock(&ctx_biglock);
  // printf("allele1: ");
  // db_nodes_print(s1, n1, db_graph, stdout);
  // printf("\nallele2: ");
  // db_nodes_print(s2, n2, db_graph, stdout);
  // printf("\n--\n");
  // pthread_mutex_unlock(&ctx_biglock);

  if(!min_covg || MIN2(sum_covg1, sum_covg2) < min_covg)
  {
    if(sum_covg1 < sum_covg2) {
      // remove s1
      mark_node_bitarr_mt(s1, n1, remove);
    }
    else {
      // remove s2
      mark_node_bitarr_mt(s2, n2, visited);
      mark_node_bitarr_mt(s2, n2, remove);
    }
  }
}

static inline void mark_remove_bubbles(dBNodeBuffer nbuf, size_t threadid,
                                       void *arg)
{
  PopBubbles *pb = (PopBubbles*)arg;
  dBNodeBuffer *alt = &pb->alts[threadid];
  const dBGraph *db_graph = pb->db_graph;

  dBNode node0, node1, nodes0[16], nodes1[16], endnode;
  uint8_t i, j, n0, n1;

  node0 = db_node_reverse(nbuf.b[0]);
  node1 = nbuf.b[nbuf.len-1];

  n0 = get_parallel_nodes(db_graph, node0, nodes0);
  n1 = get_parallel_nodes(db_graph, node1, nodes1);

  printf("n0: %zu n1: %zu\n", (size_t)n0, (size_t)n1);

  if(!n0 || !n1) return;

  for(i = 0; i < n0; i++)
  {
    db_node_buf_reset(alt);
    db_node_buf_add(alt, db_node_reverse(nodes0[i]));
    supernode_extend(alt, 0, db_graph);

    // find end node in right hand nodes
    endnode = alt->b[alt->len-1];
    for(j = 0; j < n1; j++) {
      if(db_nodes_are_equal(endnode, nodes1[j])) {
        found_bubble(nbuf.b, nbuf.len, alt->b, alt->len,
                     pb->min_covg, pb->visited, pb->remove, db_graph);
        break;
      }
    }
  }
}

// visited, remove should each have at least db_graph->capacity bits
// and should be initialised to zeros
// remove will have bits set for all nodes that should be removed
void pop_bubbles(const dBGraph *db_graph, size_t nthreads, size_t min_covg,
                 uint8_t *visited, uint8_t *remove)
{
  size_t i;
  PopBubbles data = {.visited = visited, .remove = remove,
                     .min_covg = min_covg, .db_graph = db_graph};
  data.alts = ctx_calloc(nthreads, sizeof(dBNodeBuffer));

  for(i = 0; i < nthreads; i++) db_node_buf_alloc(&data.alts[i], 256);
  supernodes_iterate(nthreads, visited, db_graph, mark_remove_bubbles, &data);
  for(i = 0; i < nthreads; i++) db_node_buf_dealloc(&data.alts[i]);
  ctx_free(data.alts);
}
