#ifndef PATH_HASH_H_
#define PATH_HASH_H_

#include "packed_path.h"
#include "hash_table.h"

#define PATH_HASH_EMPTY {.table = NULL, .num_of_buckets = 0, .bucket_size = 0, \
                         .capacity = 0, .mask = 0, .num_entries = 0}

#define PATH_HASH_UNSET (0xffffffffff)
#define PATH_HASH_ENTRY_ASSIGNED(x) ((x).hkey != PATH_HASH_UNSET)

// Packed structure is 10 bytes
// Do not use pointes to fields in this struct - they are not aligned
struct KPEntryStruct
{
  // 5 bytes each
  hkey_t hkey:40;
  PathIndex pindex:40;
} __attribute((packed));

typedef struct KPEntryStruct KPEntry;

typedef struct
{
  KPEntry *const table;
  const size_t num_of_buckets; // needs to store maximum of 1<<32
  const uint8_t bucket_size; // max value 255
  const uint64_t capacity, mask; // num_of_buckets * bucket_size
  uint8_t *const bucket_nitems; // number of items in each bucket
  uint8_t *const bktlocks; // always cast to volatile
  size_t num_entries;
} PathHash;

void path_hash_alloc(PathHash *phash, size_t mem_in_bytes);
void path_hash_dealloc(PathHash *phash);
void path_hash_reset(PathHash *phash);

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
PathIndex path_hash_get_pindex(const PathHash *phash, size_t pos);

#define PHASH_ITERATE(phash,func,...) do {                                     \
  size_t _i;                                                                   \
  for(_i = 0; _i < (phash)->capacity; _i++)                                    \
    if(PATH_HASH_ENTRY_ASSIGNED((phash)->table[_i]))                           \
      func(_i, ##__VA_ARGS__);                                                 \
} while(0)

#endif /* PATH_HASH_H_ */
