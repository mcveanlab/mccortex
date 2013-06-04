
#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

#include <inttypes.h>

#include "binary_kmer.h"

#ifndef REHASH_LIMIT
  #define REHASH_LIMIT 16
#endif
#define UNSET_BKMER (1UL<<63)
#define IDEAL_OCCUPANCY 0.85f

typedef struct
{
  BinaryKmer *const table;
  const uint64_t num_of_buckets; // needs to store maximum of 1<<32
  const uint32_t hash_mask; // this is num_of_buckets - 1
  const uint8_t bucket_size; // max value 255
  const uint64_t capacity; // num_of_buckets * bucket_size
  uint8_t *const bucket_length; // index of the next free element in bucket
  uint8_t *const bucket_fill; // number of filled entries in a bucket
  uint64_t unique_kmers;
  uint64_t collisions[REHASH_LIMIT];
} HashTable;

typedef uint64_t hkey_t;
#define HASH_NOT_FOUND UINT64_MAX

#define HASH_ENTRY_ASSIGNED(ptr) ((ptr)[0] & UNSET_BKMER)

// Get number of bytes required for a given size
size_t hash_table_mem(size_t req_capacity, size_t *act_capacity);

// Number of hash table entries for a given required capacity
size_t hash_table_cap(size_t req_capacity,
                      uint64_t *num_bckts_ptr, uint8_t *bckt_size_ptr);

// Returns NULL if not enough memory
HashTable* hash_table_alloc(HashTable *htable, uint64_t capacity);

void hash_table_dealloc(HashTable *hash_table);

hkey_t hash_table_find(const HashTable *const htable, const BinaryKmer const bkmer);
hkey_t hash_table_insert(HashTable *const htable, const BinaryKmer const bkmer);
hkey_t hash_table_find_or_insert(HashTable *htable, const BinaryKmer const bkmer,
                                 boolean *found);

void hash_table_delete(HashTable *const htable, hkey_t pos);

void hash_table_print_stats(const HashTable *const htable);
uint64_t hash_table_count_assigned_nodes(const HashTable *const htable);

#define HASH_TRAVERSE(ht,func, ...) do {                                       \
  BinaryKmer *htt_ptr = (ht)->table, *htt_end = htt_ptr + (ht)->capacity;      \
  for(; htt_ptr < htt_end; htt_ptr++) {                                        \
    if(HASH_ENTRY_ASSIGNED(*htt_ptr)) {                                        \
      func(htt_ptr - (ht)->table, ##__VA_ARGS__);                              \
    }                                                                          \
  }                                                                            \
} while(0)

#endif /* HASH_TABLE_H_ */
