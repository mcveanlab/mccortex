#ifndef POP_BUBBLES_H_
#define POP_BUBBLES_H_

#include "db_graph.h"

typedef struct
{
  int max_rmv_covg, max_rmv_klen, max_rmv_kdiff;
} PopBubblesPrefs;

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
                   uint8_t *visited, uint8_t *rmvbits);

#endif /* POP_BUBBLES_H_ */
