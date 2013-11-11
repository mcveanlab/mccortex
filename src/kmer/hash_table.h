
#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

#include <inttypes.h>

#include "binary_kmer.h"

#define REHASH_LIMIT 16
#define UNSET_BKMER (1UL<<63)
#define IDEAL_OCCUPANCY 0.75f
#define WARN_OCCUPANCY 0.9f
// bucket size must be <256
#define MAX_BUCKET_SIZE 32

typedef struct
{
  BinaryKmer *const table;
  const uint64_t num_of_buckets; // needs to store maximum of 1<<32
  const uint32_t hash_mask; // this is num_of_buckets - 1
  const uint8_t bucket_size; // max value 255
  const uint64_t capacity; // num_of_buckets * bucket_size
  // buckets[b][0] is the size of the bucket (can only increase)
  // buckets[b][1] is the number of filled entries in a bucket (can go up/down)
  uint8_t (*const buckets)[2];
  uint64_t unique_kmers;
  uint64_t collisions[REHASH_LIMIT];
} HashTable;

typedef uint64_t hkey_t;
#define HASH_NOT_FOUND UINT64_MAX

#define HASH_ENTRY_ASSIGNED(ptr) (!((ptr).b[0] & UNSET_BKMER))

// Number of hash table entries for a given required capacity
size_t hash_table_cap(size_t nkmers, boolean above_nkmers,
                      uint64_t *num_bckts_ptr, uint8_t *bckt_size_ptr);

// Get number of bytes required for a given number of kmers
// do not excess max capacity
size_t hash_table_mem(size_t nkmers, boolean above_nkmers, size_t *act_capacty_kmers);

// Returns NULL if not enough memory
HashTable* hash_table_alloc(HashTable *htable, uint64_t capacity);

void hash_table_dealloc(HashTable *hash_table);

hkey_t hash_table_find(const HashTable *const htable, const BinaryKmer bkmer);
hkey_t hash_table_insert(HashTable *const htable, const BinaryKmer bkmer);
hkey_t hash_table_find_or_insert(HashTable *htable, const BinaryKmer bkmer,
                                 boolean *found);

void hash_table_delete(HashTable *const htable, hkey_t pos);

void hash_table_print_stats(const HashTable *const htable);

// This is for debugging
uint64_t hash_table_count_assigned_nodes(const HashTable *const htable);

// Iterate over entries in the hash table
#define HASH_TRAVERSE(ht,func,...) HASH_TRAVERSE2(ht,func,##__VA_ARGS__)

// Iterate over all entries
#define HASH_TRAVERSE1(ht,func, ...) do {                                      \
  const BinaryKmer *htt_ptr = (ht)->table, *htt_end = htt_ptr + (ht)->capacity;\
  for(; htt_ptr < htt_end; htt_ptr++) {                                        \
    if(HASH_ENTRY_ASSIGNED(*htt_ptr)) {                                        \
      func(htt_ptr - (ht)->table, ##__VA_ARGS__);                              \
    }                                                                          \
  }                                                                            \
} while(0)

// Iterate over buckets, iterate over bucket contents
// Faster in low density hash tables
#define HASH_TRAVERSE2(ht,func, ...) ({                                        \
  BinaryKmer *bkt_strt = (ht)->table, *htt_ptr; size_t _b,_c;                  \
  for(_b = 0; _b < (ht)->num_of_buckets; _b++, bkt_strt += (ht)->bucket_size) {\
    for(htt_ptr = bkt_strt, _c = 0; _c < (ht)->buckets[_b][1]; htt_ptr++) {    \
      if(HASH_ENTRY_ASSIGNED(*htt_ptr)) {                                      \
        _c++; func(htt_ptr - (ht)->table, ##__VA_ARGS__);                      \
      }                                                                        \
    }                                                                          \
  }                                                                            \
})

#endif /* HASH_TABLE_H_ */
