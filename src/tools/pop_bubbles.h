#ifndef POP_BUBBLES_H_
#define POP_BUBBLES_H_

#include "db_graph.h"

// visited, rmvbits should each have at least db_graph->capacity bits
// and should be initialised to zeros
// rmvbits will have bits set for all nodes that should be removed
void pop_bubbles(const dBGraph *db_graph, size_t nthreads, size_t min_covg,
                 uint8_t *visited, uint8_t *rmvbits);

#endif /* POP_BUBBLES_H_ */
