#ifndef PATH_HASH_H_
#define PATH_HASH_H_

#include "packed_path.h"
#include "hash_table.h"

#define PATH_HASH_EMPTY {.table = NULL, .num_of_buckets = 0, .bucket_size = 0, \
                         .capacity = 0, .mask = 0, .num_entries = 0}

// typedef struct KPEntryStruct KPEntry;
typedef struct KPEntryStruct KPEntry;

struct PathHashStruct
{
  KPEntry *const table;
  const size_t num_of_buckets; // needs to store maximum of 1<<32
  const uint8_t bucket_size; // max value 255
  const uint64_t capacity, mask; // num_of_buckets * bucket_size
  size_t num_entries;
  volatile uint8_t *const bktlocks;
};

typedef struct PathHashStruct PathHash;

void path_hash_alloc(PathHash *phash, size_t mem_in_bytes);

void path_hash_dealloc(PathHash *phash);

// You must acquire the lock on the kmer before adding
// packed points to <PathLen><PackedSeq>
// *pos is set to the index of the entry if inserted or found
// Returns:
//   1  inserted
//   0  found
//  -1  out of memory
// Thread Safe: uses bucket level locks
int path_hash_find_or_insert_mt(PathHash *restrict phash, hkey_t hkey,
                                const uint8_t *restrict packed,
                                const uint8_t *restrict pstore, size_t colbytes,
                                size_t *restrict pos);

// Set pindex of newly inserted paths
void path_hash_set_pindex(PathHash *phash, size_t pos, PathIndex pindex);

// Get pindex of a path
PathIndex path_hash_get_pindex(PathHash *phash, size_t pos);

#endif /* PATH_HASH_H_ */
