#include "global.h"
#include "path_store_thread_safe.h"
#include "db_node.h"

//
// Thread safe wrapper for path_store
// path_store_mt_*() [mt for multithreaded]
//

// Returns true if added, false otherwise
// packed points to <PathLen><PackedSeq>
boolean path_store_mt_find_or_add(hkey_t kmer, dBGraph *db_graph, Colour colour,
                                  const uint8_t *packed, size_t plen)
{
  PathStore *pstore = &db_graph->pdata;
  volatile uint8_t *kmerlocks = (volatile uint8_t *)db_graph->path_kmer_locks;

// path_nbytes is the number of bytes in <PackedSeq>
  size_t path_nbytes = (plen+3)/4;

  // 1) Get lock for kmer
  // calling bitlock_yield_acquire instead of bitlock_acquire causes
  bitlock_yield_acquire(kmerlocks, kmer);

  const PathIndex next = *(volatile PathIndex*)&db_node_paths(db_graph, kmer);

  // 2) Search for path
  PathIndex match = path_store_find(pstore, next, packed, path_nbytes);
  if(match != PATH_NULL)
  {
    // => if already exist -> add colour -> release lock
    volatile uint8_t *colarr = packedpath_get_colset(pstore->store+match);
    boolean added = !bitset_get(colarr, colour);
    bitset_set(colarr, colour);
    bitlock_release(kmerlocks, kmer);
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
  uint8_t *colset = new_path+sizeof(PathIndex);

  // Prev
  packedpath_set_prev(new_path, next);
  // bitset
  memset(colset, 0, pstore->colset_bytes);
  bitset_set(colset, colour);
  // Length + Path
  memcpy(colset+pstore->colset_bytes, packed, sizeof(PathLen) + path_nbytes);

  __sync_synchronize();

  // 5) update kmer pointer
  db_node_paths(db_graph, kmer) = (uint64_t)(new_path - pstore->store);

  // Update number of kmers with paths if this the first path for this kmer
  if(next == PATH_NULL)
    __sync_add_and_fetch((volatile size_t*)&pstore->num_kmers_with_paths, 1);

  // Update number of paths
  __sync_add_and_fetch((volatile size_t*)&pstore->num_of_paths, 1);

  __sync_synchronize();

  // 6) release kmer lock
  bitlock_release(kmerlocks, kmer);

  return true;
}
