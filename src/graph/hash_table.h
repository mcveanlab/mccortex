#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

#include <inttypes.h>

#include "hash_mem.h"
#include "binary_kmer.h"

#define UNSET_BKMER_WORD (1UL<<63)

#define HT_BSIZE 0
#define HT_BITEMS 1

#define HASH_NOT_FOUND (UINT64_MAX>>1)
#define HASH_ENTRY_ASSIGNED(bkmer) (!((bkmer).b[0] & UNSET_BKMER_WORD))

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

// Safe to call on different entries at the same time
// NOT safe to do find() whilst doing delete()
void hash_table_delete(HashTable *const htable, hkey_t pos);

// Delete all entries from a hash table
void hash_table_empty(HashTable *const htable);

void hash_table_print_stats(const HashTable *const htable);
void hash_table_print_stats_brief(const HashTable *const htable);

// This is for debugging
uint64_t hash_table_count_kmers(const HashTable *const htable);

// Iterate over entries in the hash table
#define HASH_ITERATE(ht,func,...) HASH_ITERATE2(ht,func,##__VA_ARGS__)

// This iterator allows adding/removing items
#define HASH_ITERATE_SAFE(ht,func,...) HASH_ITERATE1(ht,func,##__VA_ARGS__)

// Iterate over all entries
#define HASH_ITERATE1(ht,func, ...) do {                                       \
  const BinaryKmer *htt_ptr = (ht)->table, *htt_end = htt_ptr + (ht)->capacity;\
  for(; htt_ptr < htt_end; htt_ptr++) {                                        \
    if(HASH_ENTRY_ASSIGNED(*htt_ptr)) {                                        \
      func((hkey_t)(htt_ptr - (ht)->table), ##__VA_ARGS__);                    \
    }                                                                          \
  }                                                                            \
} while(0)

// Iterate over buckets, iterate over bucket contents
// Faster in low density hash tables
// Don't use this iterator if your func adds or removes elements
#define HASH_ITERATE2(ht,func, ...) do {                                       \
  const BinaryKmer *bkt_strt = (ht)->table, *htt_ptr; size_t _b,_c;            \
  for(_b = 0; _b < (ht)->num_of_buckets; _b++, bkt_strt += (ht)->bucket_size) {\
    for(htt_ptr = bkt_strt, _c = 0; _c < (ht)->buckets[_b][HT_BITEMS]; htt_ptr++){\
      if(HASH_ENTRY_ASSIGNED(*htt_ptr)) {                                      \
        _c++; func((hkey_t)(htt_ptr - (ht)->table), ##__VA_ARGS__);            \
      }                                                                        \
    }                                                                          \
  }                                                                            \
} while(0)

// This iterator allows adding/removing items
#define HASH_ITERATE_PART(ht,job,njobs,func, ...) do {                         \
  const size_t _step = (ht)->capacity / (njobs);                               \
  const BinaryKmer *_start, *_end, *_bkptr;                                    \
  _start = (ht)->table + (job) * _step;                                        \
  _end = ((job)+1 == (njobs) ? (ht)->table + (ht)->capacity : _start+_step);   \
  for(_bkptr = _start; _bkptr < _end; _bkptr++) {                              \
    if(HASH_ENTRY_ASSIGNED(*_bkptr)) {                                         \
      func((hkey_t)(_bkptr - (ht)->table), ##__VA_ARGS__);                     \
    }                                                                          \
  }                                                                            \
} while(0)

#endif /* HASH_TABLE_H_ */
