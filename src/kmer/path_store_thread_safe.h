#ifndef PATH_STORE_THREAD_SAFE_H_
#define PATH_STORE_THREAD_SAFE_H_

#include "db_graph.h"
#include "path_store.h"

// Returns true if new to colour, false otherwise
// packed points to <PathLen><PackedSeq>
// Returns address of path in PathStore by setting newidx
boolean path_store_mt_find_or_add(hkey_t hkey, dBGraph *db_graph, Colour ctpcol,
                                  const uint8_t *packed, size_t plen,
                                  PathIndex *newidx);

#endif /* PATH_STORE_THREAD_SAFE_H_ */
