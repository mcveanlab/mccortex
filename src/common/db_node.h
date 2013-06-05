#ifndef DB_NODE_H_
#define DB_NODE_H_

#include <inttypes.h>

#include "binary_kmer.h"
#include "db_graph.h"

#define COVG_MAX UINT_MAX

//
// Status Constants
//
#define EFLAG_ZERO      0x0

// Removing duplicates when loading sequence (forward and reverse)
#define EFLAG_RDSTRT_FW 0x4 // bit:2
#define EFLAG_RDSTRT_RV 0x8 // bit:3
// Walking the graph
#define EFLAG_WALK_FW   0x10// bit:4
#define EFLAG_WALK_RV   0x20// bit:5
// Dump a subgraph
#define EFLAG_MARKED    0x40// bit:6
// Adding shades
#define EFLAG_IN_READ   0x40// bit:6

//
// Status
//
#define db_node_set_flag(graph,hkey,flag)  ((graph)->status[hkey] |= (flag))
#define db_node_del_flag(graph,hkey,flag)  ((graph)->status[hkey] &=~(flag))
#define db_node_has_flag(graph,hkey,flag)  ((graph)->status[hkey] &  (flag))
#define db_node_reset_flags(graph,hkey)    ((graph)->status[hkey] = 0)

//
// Get binary kmer key
//
Key db_node_get_key(const uint64_t* const restrict kmer, uint32_t kmer_size,
                    Key restrict kmer_key);

#define db_node_to_str(graph,node,str) \
        binary_kmer_to_str(db_graph_bkmer(graph,node), (graph)->kmer_size, (str))

//
// Orientations
//
#define rev_orient(or) (!(or))
#define opposite_orientation(or) rev_orient(or)

#define db_node_get_orientation(bkmer,bkey) \
        (binary_kmers_are_equal((bkmer), (bkey)) ? forward : reverse)

void db_node_oriented_bkmer(const BinaryKmer bkmer, Orientation orient,
                            uint32_t kmer_size, BinaryKmer result);

#define db_node_first_nuc(bkmer,or,k) \
  ((or) == forward ? binary_kmer_first_nuc((bkmer),(k)) \
      : binary_nuc_complement(binary_kmer_last_nuc(bkmer)))

#define db_node_last_nuc(bkmer,or,k) \
  ((or) == forward ? binary_kmer_last_nuc(bkmer) \
      : binary_nuc_complement(binary_kmer_first_nuc(bkmer,(k))))

//
// Edges
//

// Orientation is 0(forward) or 1(reverse), shifted left (<<) 2 gives 0 or 4
#define nuc_orient_to_edge(n,or)    (0x1 << ((n) + ((or)<<2)))

#define edges_has_edge(edges,n,or)  (((edges) >> ((n) + ((or)<<2))) & 0x1)
#define edges_set_edge(edges,n,or)  ((edges) | nuc_orient_to_edge(n,or))
#define edges_del_edge(edges,n,or)  ((edges) &~nuc_orient_to_edge(n,or))

#define edges_with_orientation(edges,or) (((edges) >> ((or)<<2)) & 0xf)

#define edges_get_outdegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,or))

#define edges_get_indegree(edges,or) \
        __builtin_popcount(edges_with_orientation(edges,rev_orient(or)))


#define db_node_has_edge(graph,hkey,nuc,or) \
        edges_has_edge((graph)->edges[hkey],(nuc),(or))
#define db_node_del_edge(graph,hkey,nuc,or) \
        ((graph)->edges[hkey] = edges_del_edge((graph)->edges[hkey],(nuc),(or)))
#define db_node_set_edge(graph,hkey,nuc,or) \
        ((graph)->edges[hkey] = edges_set_edge((graph)->edges[hkey],(nuc),(or)))

#define db_node_get_edges(graph,hkey) ((graph)->edges[hkey])
#define db_node_reset_edges(graph,hkey) ((graph)->edges[hkey] = 0)

#define db_node_is_blunt_end(graph,hkey,or) \
        (edges_with_orientation((graph)->edges[hkey],or) == 0)

boolean edges_has_precisely_one_edge(Edges edges, Orientation orientation,
                                     Nucleotide *nucleotide);

//
// Coverages
//

#define db_node_set_covg(graph,hkey,col,covg) ((graph)->covgs[col][hkey] = (covg))
#define db_node_get_covg(graph,hkey,col)      ((graph)->covgs[col][hkey])

void db_node_add_coverage(dBGraph *graph, hkey_t hkey, Colour col, long update);
void db_node_increment_coverage(dBGraph *graph, hkey_t hkey, Colour col);

Covg db_node_sum_covg_of_colours(const dBGraph *graph, hkey_t hkey,
                                 Colour first_col, int num_cols);
Covg db_node_sum_covg_of_all_colours(const dBGraph *graph, hkey_t hkey);
Covg db_node_sum_covg_of_colourlist(const dBGraph *graph, hkey_t hkey,
                                    const Colour *colour_list, int num_cols);

//
// dBNodeBuffer
//
typedef struct
{
  hkey_t *nodes;
  Orientation *orients;
  size_t len, capacity;
} dBNodeBuffer;

void db_node_buf_alloc(dBNodeBuffer *buf, size_t capacity);
void db_node_buf_dealloc(dBNodeBuffer *buf);
void db_node_buf_free(dBNodeBuffer *buf);
void db_node_buf_ensure_capacity(dBNodeBuffer *buf, size_t capacity);
void db_node_buf_safe_add(dBNodeBuffer *buf, hkey_t node, Orientation orient);

#endif /* DB_NODE_H_ */
