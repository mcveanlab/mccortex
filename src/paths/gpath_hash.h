#ifndef GPATH_HASH_H_
#define GPATH_HASH_H_

#include "cortex_types.h"
#include "gpath_store.h"

// Packed structure is 10 bytes
// Do not use pointes to fields in this struct - they are not aligned
struct GPEntryStruct
{
  // 5 bytes each
  hkey_t hkey:40;
  pkey_t gpindex:40;
} __attribute((packed));

/*
// 5+5+5+1+2 = 18 bytes
// GPath:18 + GPEntry:10 + count:1 + colset:1 = 30
// 18/34 = 52%
struct GPEntryStruct
{
  hkey_t hkey:40; // 5 bytes
  pkey_t next:40; // 5 bytes
  uint64_t seq_offset:40; // 5 btyes
  uint8_t count; // 1 byte
  uint16_t orient:1, seq_len:15; // 2 bytes
} __attribute((packed));
*/

typedef struct GPEntryStruct GPEntry;

typedef struct
{
  GPathStore *const gpstore; // Add to this path store
  GPEntry *const table; // Using this table to remove duplicates
  const size_t num_of_buckets; // needs to store maximum of 1<<32
  const uint8_t bucket_size; // max value 255
  const uint64_t capacity, mask; // num_of_buckets * bucket_size
  uint8_t *const bucket_nitems; // number of items in each bucket
  uint8_t *const bktlocks; // always cast to volatile
  // uint8_t *const seq;
  // size_t seq_len, seq_capacity;
  size_t num_entries;
} GPathHash;

void gpath_hash_alloc(GPathHash *phash, GPathStore *gpstore, size_t mem_in_bytes);
void gpath_hash_dealloc(GPathHash *phash);
void gpath_hash_reset(GPathHash *phash);

void gpath_hash_print_stats(const GPathHash *phash);

// Returns NULL if out of memory
// Thread Safe: uses bucket level locks
GPath* gpath_hash_find_or_insert_mt(GPathHash *restrict phash,
                                    hkey_t hkey, GPathNew newgpath,
                                    bool *found);

#endif /* GPATH_HASH_H_ */
