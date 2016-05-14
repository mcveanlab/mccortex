#include "global.h"
#include "pop_bubbles.h"
#include "db_unitig.h"

/*
  Popping bubbles works by iterating over all unitigs. For each unitig
  we attempt to pull out parallel unitigs. Once we have two parallel
  unitigs we see if one should be removed.
*/

/*
 In some rare situations boths branches of a bubble will be removed.
 Example: Both bottom branches may be removed.

   ___         ___         ___         ___  
 _/   \_     _/   \_     _/ ^ \_     _/ ^ \_
  \___/       \___/       \___/             
 _/ ^ \_ ->  _/   \_ ->  _/ ^ \_ ->  _/   \_
  \___/                                     
    ^                                       

  Note: Could avoid this if we only remove unitigs
        where in-degree == 1 and out-degree == 1,
        otherwise just remove our edges.

   ___         ___         ___         ___  
 _/   \_     _/   \_     _/ ^ \_     _/ ^ \_
  \___/       \___/       \___/        ___  
 _/ ^ \_ ->  _/   \_ ->  _/ ^ \_ ->  _/   \_
  \___/                                     
    ^                                       

 */

typedef struct
{
  uint8_t *const visited, *const rmvbits;
  dBNodeBuffer *alts;
  size_t *num_popped;
  const PopBubblesPrefs prefs;
  const dBGraph *db_graph;
} PopBubbles;

/*
  Go forward one node and back one node to get all 'parallel nodes'
  Given node, return nodes {a,b}
       a
        \
  node -> x
        /
       b
*/
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
  for(i = 0; i < n; i++) (void)bitset_set_mt(bitset, nodes[i].key);
}

/**
 * Remove the lowest mean coverage branch, by marking the remove bit array.
 * @param min_covg keep all branches with mean coverage >= min_covg
 */
static inline bool process_bubble(const dBNode *s1, size_t n1,
                                  const dBNode *s2, size_t n2,
                                  const PopBubblesPrefs *p,
                                  uint8_t *visited, uint8_t *rmvbits,
                                  const dBGraph *db_graph)
{
  size_t i, sum_covg1 = 0, sum_covg2 = 0, mean_covg1, mean_covg2;

  for(i = 0; i < n1; i++) sum_covg1 += db_node_sum_covg(db_graph, s1[i].key);
  for(i = 0; i < n2; i++) sum_covg2 += db_node_sum_covg(db_graph, s2[i].key);
  mean_covg1 = sum_covg1 / n1;
  mean_covg2 = sum_covg2 / n2;

  // Print alleles
  // pthread_mutex_lock(&ctx_biglock);
  // printf("allele1: ");
  // db_nodes_print(s1, n1, db_graph, stdout);
  // printf("\nallele2: ");
  // db_nodes_print(s2, n2, db_graph, stdout);
  // printf("\n--\n");
  // pthread_mutex_unlock(&ctx_biglock);

  size_t rmv_covg, rmv_klen;
  if(mean_covg1 < mean_covg2) { rmv_covg = mean_covg1; rmv_klen = n1; }
  else                        { rmv_covg = mean_covg2; rmv_klen = n2; }

  if((!p->max_rmv_covg     || rmv_covg <= (size_t)p->max_rmv_covg) &&
     (!p->max_rmv_klen     || rmv_klen <= (size_t)p->max_rmv_klen) &&
     (p->max_rmv_kdiff < 0 || abs((int)n1 - (int)n2) <= p->max_rmv_kdiff))
  {
    if(mean_covg1 < mean_covg2) {
      // remove s1
      mark_node_bitarr_mt(s1, n1, rmvbits);
    }
    else {
      // remove s2
      mark_node_bitarr_mt(s2, n2, visited);
      mark_node_bitarr_mt(s2, n2, rmvbits);
    }
    return true;
  }
  return false;
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

  // printf("n0: %zu n1: %zu\n", (size_t)n0, (size_t)n1);

  if(!n0 || !n1) return;

  for(i = 0; i < n0; i++)
  {
    db_node_buf_reset(alt);
    db_node_buf_add(alt, db_node_reverse(nodes0[i]));
    db_unitig_extend(alt, 0, db_graph);

    // find end node in right hand nodes
    endnode = alt->b[alt->len-1];
    for(j = 0; j < n1; j++) {
      if(db_nodes_are_equal(endnode, nodes1[j])) {
        // found a bubble
        if(process_bubble(nbuf.b, nbuf.len, alt->b, alt->len,
                          &pb->prefs, pb->visited, pb->rmvbits, db_graph))
        {
          // Popped a bubble
          pb->num_popped[threadid]++;
        }
        break;
      }
    }
  }
}

/**
 * visited, rmvbits should each have at least db_graph->capacity bits
 * and should be initialised to zeros
 * rmvbits will have bits set for all nodes that should be removed
 * @param max_rmv_covg only remove contigs with covg <= max_rmv_covg,
 *                     ignored if <= 0.
 * @param max_rmv_klen only remove contigs with num kmers <= max_rmv_klen,
 *                     ignored if <= 0.
 * @param max_rmv_kdiff only remove contigs if max diff in kmers <= max_rmv_kdiff,
 *                      ignored if < 0.
 * @return number of bubbles popped
**/
size_t pop_bubbles(const dBGraph *db_graph, size_t nthreads,
                   PopBubblesPrefs prefs,
                   uint8_t *visited, uint8_t *rmvbits)
{
  size_t i, total_popped = 0;

  status("[pop_bubbles] Popping bubbles...");
  if(prefs.max_rmv_covg > 0)
    status("[pop_bubbles]   where branch coverage <= %i", prefs.max_rmv_covg);
  if(prefs.max_rmv_klen > 0)
    status("[pop_bubbles]   where branch length <= %i", prefs.max_rmv_klen);
  if(prefs.max_rmv_kdiff >= 0)
    status("[pop_bubbles]   where branch length diff < %i", prefs.max_rmv_kdiff);

  PopBubbles data = {.visited = visited, .rmvbits = rmvbits,
                     .prefs = prefs, .db_graph = db_graph};

  data.alts = ctx_calloc(nthreads, sizeof(dBNodeBuffer));
  data.num_popped = ctx_calloc(nthreads, sizeof(size_t));
  for(i = 0; i < nthreads; i++) db_node_buf_alloc(&data.alts[i], 256);

  db_unitigs_iterate(nthreads, visited, db_graph, mark_remove_bubbles, &data);

  for(i = 0; i < nthreads; i++) {
    total_popped += data.num_popped[i];
    db_node_buf_dealloc(&data.alts[i]);
  }
  ctx_free(data.num_popped);
  ctx_free(data.alts);

  return total_popped;
}
