#include "global.h"
#include "path_store_thread_safe.h"
#include "db_node.h"
#include "graph_paths.h"

//
// Thread safe wrapper for path_store
// path_store_mt_*() [mt for multithreaded]
//

// Returns true if new to colour, false otherwise
// packed points to <PathLen><PackedSeq>
// Returns address of path in PathStore by setting newidx
boolean path_store_mt_find_or_add(hkey_t hkey, dBGraph *db_graph, Colour ctpcol,
                                  const uint8_t *packed, size_t plen,
                                  PathIndex *newidx)
{
  PathStore *pstore = &db_graph->pdata;
  volatile uint8_t *kmerlocks = (volatile uint8_t *)db_graph->path_kmer_locks;

// path_nbytes is the number of bytes in <PackedSeq>
  size_t path_nbytes = (plen+3)/4;

  // debug
  // print_path(hkey, packed, pstore);

  // 1) Get lock for kmer
  // calling bitlock_yield_acquire instead of bitlock_acquire causes
  bitlock_yield_acquire(kmerlocks, hkey);

  const PathIndex next = *(volatile PathIndex*)&db_node_paths(db_graph, hkey);

  // 2) Search for path
  PathIndex match = path_store_find(pstore, next, packed, path_nbytes);
  uint8_t *colset;

  if(match != PATH_NULL)
  {
    // => if already exist -> add colour -> release lock
    colset = packedpath_get_colset(pstore->store+match);
    boolean added = !bitset_get(colset, ctpcol);
    bitset_set(colset, ctpcol);
    bitlock_release(kmerlocks, hkey);
    *newidx = match;
    return added;
  }

  // 3) shift store->next to make space
  size_t mem = packedpath_mem2(pstore->colset_bytes, path_nbytes);
  uint8_t *new_path;

  // atomic { new_path = pstore->next; pstore->next += mem; }
  // __sync_fetch_and_add(x,y) x and y need to be of same type
  new_path = __sync_fetch_and_add((uint8_t*volatile*)&pstore->next,
                                  (uint8_t*)mem);

  if(new_path + mem > pstore->end) die("Out of path memory!");

  // 4) Copy new entry

  // Prev
  packedpath_set_prev(new_path, next);

  // bitset
  colset = packedpath_get_colset(new_path);
  memset(colset, 0, pstore->colset_bytes);
  bitset_set(colset, ctpcol);

  // Length + Path
  memcpy(colset+pstore->colset_bytes, packed, sizeof(PathLen) + path_nbytes);

  // path must be written before we move the kmer path pointer forward
  // although there is a write-lock (kmerlocks), threads currently traversing
  // the graph would fall over
  __sync_synchronize();

  // 5) update kmer pointer
  PathIndex pindex = (uint64_t)(new_path - pstore->store);
  db_node_paths(db_graph, hkey) = pindex;

  // Update number of kmers with paths if this the first path for this kmer
  if(next == PATH_NULL)
    __sync_add_and_fetch((volatile size_t*)&pstore->num_kmers_with_paths, 1);

  // Update number of paths
  __sync_add_and_fetch((volatile size_t*)&pstore->num_of_paths, 1);

  __sync_synchronize();

  // 6) release kmer lock
  bitlock_release(kmerlocks, hkey);

  *newidx = pindex;
  return true;
}
