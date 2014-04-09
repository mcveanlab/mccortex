#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

#include <inttypes.h>

#include "binary_kmer.h"

#define REHASH_LIMIT 20
#define UNSET_BKMER_WORD (1UL<<63)
#define IDEAL_OCCUPANCY 0.75f
#define WARN_OCCUPANCY 0.9f
// bucket size must be <256
#define MAX_BUCKET_SIZE 48

#define HT_BSIZE 0
#define HT_BITEMS 1

// Struct is public so ITERATE macros can operate on it
typedef struct
{
  BinaryKmer *const table;
  const uint64_t num_of_buckets; // needs to store maximum of 1<<32
  const uint_fast32_t hash_mask; // this is num_of_buckets - 1
  const uint8_t bucket_size; // max value 255
  const uint64_t capacity; // num_of_buckets * bucket_size
  // buckets[b][0] is the size of the bucket (can only increase)
  // buckets[b][1] is the number of filled entries in a bucket (can go up/down)
  uint8_t (*const buckets)[2];
  uint64_t num_kmers;
  uint64_t collisions[REHASH_LIMIT];
} HashTable;

typedef uint64_t hkey_t; // don't ever use the top bit, used later for orientation

typedef struct {
  hkey_t orient:1, key:63;
} dBNode;

#define HASH_NOT_FOUND (UINT64_MAX>>1)

#define HASH_ENTRY_ASSIGNED(ptr) (!((ptr).b[0] & UNSET_BKMER_WORD))

// Hash table capacity is x*(2^y) where x and y are parameters
// memory is x*(2^y)*sizeof(BinaryKmer) + (2^y) * 2
#define ht_mem(bktsize,nbkts,nbits) \
        (((bktsize) * (nbkts) * (sizeof(BinaryKmer)*8+(nbits)))/8 +\
         (nbkts) * sizeof(uint8_t[2]))

// Returns capacity of a hash table that holds at least nkmers
size_t hash_table_cap(size_t nkmers, uint64_t *num_bkts_ptr, uint8_t *bkt_size_ptr);

// Returns memory required to hold nkmers
size_t hash_table_mem(size_t nkmers, size_t extrabits, size_t *nkmers_ptr);

// Returns memory used for hashtable no more than some memory limit
size_t hash_table_mem_limit(size_t memlimit, size_t extrabits, size_t *nkmers_ptr);

// Returns NULL if not enough memory
void hash_table_alloc(HashTable *htable, uint64_t capacity);
void hash_table_dealloc(HashTable *hash_table);

hkey_t hash_table_find(const HashTable *const htable, const BinaryKmer bkmer);
hkey_t hash_table_insert(HashTable *const htable, const BinaryKmer bkmer);
hkey_t hash_table_find_or_insert(HashTable *htable, const BinaryKmer bkmer,
                                 bool *found);

// Threadsafe find or insert, using bucket level locks
hkey_t hash_table_find_or_insert_mt(HashTable *htable, const BinaryKmer key,
                                    bool *found, volatile uint8_t *bktlocks);

void hash_table_delete(HashTable *const htable, hkey_t pos);

// Delete all entries from a hash table
void hash_table_empty(HashTable *const htable);

void hash_table_print_stats(const HashTable *const htable);
void hash_table_print_stats_brief(const HashTable *const htable);

// This is for debugging
uint64_t hash_table_count_kmers(const HashTable *const htable);

// Iterate over entries in the hash table
#define HASH_ITERATE(ht,func,...) HASH_ITERATE2(ht,func,##__VA_ARGS__)

// This allows up to add/remove items
#define HASH_ITERATE_SAFE(ht,func,...) HASH_ITERATE1(ht,func,##__VA_ARGS__)

// Iterate over all entries
#define HASH_ITERATE1(ht,func, ...) {                                          \
  const BinaryKmer *htt_ptr = (ht)->table, *htt_end = htt_ptr + (ht)->capacity;\
  for(; htt_ptr < htt_end; htt_ptr++) {                                        \
    if(HASH_ENTRY_ASSIGNED(*htt_ptr)) {                                        \
      func((hkey_t)(htt_ptr - (ht)->table), ##__VA_ARGS__);                    \
    }                                                                          \
  }                                                                            \
}

// Iterate over buckets, iterate over bucket contents
// Faster in low density hash tables
// Don't use this iterator if your func adds or removes elements
#define HASH_ITERATE2(ht,func, ...) {                                          \
  const BinaryKmer *bkt_strt = (ht)->table, *htt_ptr; size_t _b,_c;            \
  for(_b = 0; _b < (ht)->num_of_buckets; _b++, bkt_strt += (ht)->bucket_size) {\
    for(htt_ptr = bkt_strt, _c = 0; _c < (ht)->buckets[_b][HT_BITEMS]; htt_ptr++){\
      if(HASH_ENTRY_ASSIGNED(*htt_ptr)) {                                      \
        _c++; func((hkey_t)(htt_ptr - (ht)->table), ##__VA_ARGS__);            \
      }                                                                        \
    }                                                                          \
  }                                                                            \
}

#endif /* HASH_TABLE_H_ */
