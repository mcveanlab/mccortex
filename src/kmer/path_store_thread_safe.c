#include "global.h"
#include "path_store_thread_safe.h"
#include "db_node.h"

//
// Thread safe wrapper for path_store
// path_store_mt_*() [mt for multithreaded]
//

// Returns true if added, false otherwise
boolean path_store_mt_find_or_add(hkey_t kmer, dBGraph *db_graph, Colour colour,
                                  const uint8_t *packed, size_t path_nbytes)
{
  PathStore *pstore = &db_graph->pdata;
  volatile uint8_t *kmerlocks = (volatile uint8_t *)db_graph->path_kmer_locks;

  // 1) Get lock for kmer
  bitlock_acquire(kmerlocks, kmer);

  PathIndex next = *(volatile PathIndex*)&db_node_paths(db_graph, kmer);

  // 2) Search for path
  PathIndex match = path_store_find(pstore, next, packed, path_nbytes);
  if(match != PATH_NULL)
  {
    // => if already exist -> add colour -> release lock
    volatile uint8_t *colarr = packedpath_colset(pstore->store+match);
    boolean added = !bitset_get(colarr, colour);
    bitset_set(colarr, colour);
    bitlock_release(kmerlocks, kmer);
    return added;
  }

  // 3) shift store->next to make space
  size_t mem = packedpath_mem2(pstore->colset_bytes, path_nbytes);
  uint8_t *write;
  // write = pstore->next; pstore->next += mem;
  // __sync_fetch_and_add(x,y) x and y need to be of same type
  write = (uint8_t*)__sync_fetch_and_add((volatile uint8_t**)&pstore->next,
                                         (uint8_t*)mem);

  if(write + mem > pstore->end) die("Out of path memory!");

  // 4) Copy new entry

  // Prev
  packedpath_prev(write) = next;
  write += sizeof(PathIndex);
  // bitset
  memset(write, 0, pstore->colset_bytes);
  bitset_set(write, colour);
  write += pstore->colset_bytes;
  // Length
  *(PathLen*)write = (PathLen)path_nbytes;
  write += sizeof(PathLen);
  // Path
  memcpy(write, packed, path_nbytes);

  __sync_synchronize();

  // 5) update kmer pointer
  db_node_paths(db_graph, kmer) = (uint64_t)(write - pstore->store);

  __sync_synchronize();

  // 6) release kmer lock
  bitlock_release(kmerlocks, kmer);

  return true;
}
