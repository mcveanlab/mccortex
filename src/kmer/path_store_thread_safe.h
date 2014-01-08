#ifndef PATH_STORE_THREAD_SAFE_H_
#define PATH_STORE_THREAD_SAFE_H_

#include "db_graph.h"
#include "path_store.h"

boolean path_store_mt_find_or_add(hkey_t kmer, dBGraph *db_graph, Colour colour,
                                  const uint8_t *packed, size_t path_nbytes);

#endif /* PATH_STORE_THREAD_SAFE_H_ */
