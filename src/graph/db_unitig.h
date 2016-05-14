#ifndef UNITIG_H_
#define UNITIG_H_

#include "cortex_types.h"
#include "db_graph.h"
#include "db_node.h"

// Orient unitig
// Once oriented, unitig has lowest poosible kmerkey at the beginning,
// oriented FORWARDs if possible
void db_unitig_normalise(dBNode *nlist, size_t len, const dBGraph *db_graph);

// Extend a unitig, nlist[offset] and olist[offset] must already be set
// Walk along nodes starting from node/or, storing the unitig in nlist/olist
// Returns the number of nodes added, adds no more than `limit`
// return false if out of space and limit > 0
bool db_unitig_extend(dBNodeBuffer *nbuf, size_t limit,
                      const dBGraph *db_graph);

// Fills with unitig that contains hkey
// Does not reset nbuf
void db_unitig_fetch(hkey_t node, dBNodeBuffer *nbuf, const dBGraph *db_graph);

// Count number of read starts using coverage data
size_t db_unitig_read_starts(const Covg *covgs, size_t len);
size_t db_unitig_covg_mean(const Covg *covgs, size_t len);

/**
 * @param visited must be initialised to zero, will be dirty upon return
 **/
void db_unitigs_iterate(size_t nthreads, uint8_t *visited,
                        const dBGraph *db_graph,
                        void (*func)(dBNodeBuffer nbuf, size_t threadid, void *arg),
                        void *arg);

#endif /* UNITIG_H_ */
