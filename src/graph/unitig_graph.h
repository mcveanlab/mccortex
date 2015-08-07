#ifndef UNITIG_GRAPH_H_
#define UNITIG_GRAPH_H_

#include "db_graph.h"

/**
 * Each unitig end is packed into 64 bits
 * Each kmer on the end of a unitig has unitig info
 * @field assigned is true iff this kmer is the first and/or last kmer of a unitig
 * @field unitigid is the index of the unitig the given kmer belongs to
 * @field left is true iff given kmer is the first kmer in a unitig
 * @field right is true iff given kmer is the last kmer in a unitig
 * @field lorient is the orientation of the first kmer at the start of the unitig
 * @field rorient is the orientation of the last kmer at the end of the unitig
 */
typedef struct {
  size_t unitigid:57, left:1, right:1, lorient:1, rorient:1, assigned:1;
} UnitigEnd;

// A UnitigKmerGraph labels each kmer at the end of a unitig
typedef struct {
  UnitigEnd *unitig_ends;
  size_t num_unitigs;
  const dBGraph *db_graph;

  // If set, during construction function is called on each unitig
  void (*per_untig)(const dBNode *nodes, size_t n, size_t uidx, void *arg);
  void *per_untig_arg;
} UnitigKmerGraph;

// Returns unitig ID
// threadsafe
size_t unitig_graph_store_end_mt(const dBNode *nodes, size_t n,
                                 UnitigKmerGraph *ugraph);

/**
 * @param visited must be initialised to zero, will be dirty upon return
 **/
void unitig_graph_create(UnitigKmerGraph *ugraph,
                         size_t nthreads,
                         uint8_t *visited,
                         void (*per_untig)(const dBNode *nodes, size_t n,
                                           size_t uidx, void *arg),
                        void *per_untig_arg);

void unitig_graph_alloc(UnitigKmerGraph *ugraph, const dBGraph *db_graph);
void unitig_graph_dealloc(UnitigKmerGraph *ugraph);

#endif /* UNITIG_GRAPH_H_ */
