#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

#include <inttypes.h>

#include "hash_mem.h"
#include "binary_kmer.h"
#include "util.h"

#define HT_BSIZE 0
#define HT_BITEMS 1

#define HASH_NOT_FOUND (UINT64_MAX>>1)
#define BKMER_SET_FLAG (1UL<<63)
#define HASH_ENTRY_ASSIGNED(bkmer) (((bkmer).b[0] & BKMER_SET_FLAG))

// Struct is public so ITERATE macros can operate on it
typedef struct
{
  BinaryKmer *const table; // Do not directly access use hash_table_fetch()!
  const uint64_t num_of_buckets; // needs to store maximum of 1<<32
  const uint_fast32_t hash_mask; // this is num_of_buckets - 1
  const uint8_t bucket_size; // max value 255
  const uint64_t capacity; // num_of_buckets * bucket_size
  // buckets[b][0] is the size of the bucket (can only increase)
  // buckets[b][1] is the number of filled entries in a bucket (can go up/down)
  uint8_t (*const buckets)[2];
  uint64_t num_kmers;
  uint64_t collisions[REHASH_LIMIT];
  const uint32_t seed; // random seed used in hashing
} HashTable;

// Returns NULL if not enough memory
void hash_table_alloc(HashTable *htable, uint64_t capacity);
void hash_table_dealloc(HashTable *ht);

#define hash_table_size(ht) (ht)->capacity
#define hash_table_nkmers(ht) (ht)->num_kmers
#define hash_table_assigned(ht,key) HASH_ENTRY_ASSIGNED((ht)->table[key])

static inline BinaryKmer hash_table_fetch(const HashTable *const ht, hkey_t key)
{
  BinaryKmer bk = ht->table[key];
  bk.b[0] &= 0x3fffffffffffffff; // mask off top two bits
  return bk;
}

#define hash_table_nbuckets(ht) ((ht)->num_of_buckets)
#define hash_table_bucket_size(ht) ((ht)->bucket_size)
#define hash_table_bsize(ht,bkt) ((ht)->buckets[bkt][HT_BSIZE])
#define hash_table_bitems(ht,bkt) ((ht)->buckets[bkt][HT_BITEMS])

hkey_t hash_table_find(const HashTable *const htable, const BinaryKmer bkmer);
hkey_t hash_table_insert(HashTable *const htable, const BinaryKmer bkmer);
hkey_t hash_table_find_or_insert(HashTable *htable, const BinaryKmer bkmer,
                                 bool *found);

// Threadsafe find, using bucket level locks
hkey_t hash_table_find_mt(HashTable *ht, const BinaryKmer key,
                          volatile uint8_t *bktlocks);

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

// Returns sorted array of hkey_t from the hash table, use kmers[i].h
hkey_t* hash_table_sorted(const HashTable *htable);

// This is for debugging
uint64_t hash_table_count_kmers(const HashTable *const htable);

// Iterate over entries in the hash table
// calls func(hkey,...)
#define HASH_ITERATE(ht,func,...) HASH_ITERATE2(ht,func,##__VA_ARGS__)

// This iterator allows adding/removing items
#define HASH_ITERATE_SAFE(ht,func,...) HASH_ITERATE1(ht,func,##__VA_ARGS__)

// Iterate over all entries
#define HASH_ITERATE1(ht,func, ...) do {                                       \
  hkey_t _hi, _hsize = hash_table_size(ht);                                    \
  for(_hi = 0; _hi < _hsize; _hi++) {                                          \
    if(hash_table_assigned(ht, _hi)) {                                         \
      func(_hi, ##__VA_ARGS__);                                                \
    }                                                                          \
  }                                                                            \
} while(0)

// Iterate over buckets, iterate over bucket contents
// Faster in low density hash tables
// Don't use this iterator if your func adds or removes elements
#define HASH_ITERATE2(ht,func, ...) do {                                       \
  size_t _b, _bstart, _nitems, _nseen, _hi;                                    \
  size_t _nbuck = hash_table_nbuckets(ht), _bsize = hash_table_bucket_size(ht);\
  for(_b = _bstart = 0; _b < _nbuck; _b++, _bstart += _bsize) {                \
    _nitems = hash_table_bitems(ht, _b);                                       \
    for(_nseen = 0, _hi = _bstart; _nseen < _nitems; _hi++) {                  \
      if(hash_table_assigned(ht, _hi)) { _nseen++; func(_hi, ##__VA_ARGS__); } \
    }                                                                          \
  }                                                                            \
} while(0)

// iterate over kmers in lexigraphic order
// Requires sizeof(hkey_t) * ht->num_kmers memory which it allocates and frees
#define HASH_ITERATE_SORTED(ht,func, ...) do {                                 \
  hkey_t *_hkeys = hash_table_sorted(ht);                                      \
  size_t _i, _nkmers = hash_table_nkmers(ht);                                  \
  for(_i = 0; _i < _nkmers; _i++) { func(_hkeys[_i], ##__VA_ARGS__); }         \
  ctx_free(_hkeys);                                                            \
} while(0)

// This iterator allows adding/removing items
// Stops if func() returns non-zero value
#define HASH_ITERATE_PART(ht,job,njobs,func, ...) do {                         \
  ctx_assert((job) < (njobs));                                                 \
  const size_t _step = hash_table_size(ht) / (njobs);                          \
  hkey_t _hi = (job) * _step;                                                  \
  hkey_t _end = ((job)+1 == (njobs) ? hash_table_size(ht) : _hi+_step);        \
  for(; _hi < _end; _hi++) {                                                   \
    if(hash_table_assigned(ht,_hi)) {                                          \
      if(func(_hi, ##__VA_ARGS__)) break;                                      \
    }                                                                          \
  }                                                                            \
} while(0)

typedef struct
{
  const HashTable *const ht;
  const size_t nthreads;
  bool (*const func)(hkey_t _h, size_t threadid, void *_arg);
  void *arg;
} HashTableIterator;

static inline void _hash_table_iterate(void *arg, size_t threadid)
{
  HashTableIterator itr = *(HashTableIterator*)arg;
  HASH_ITERATE_PART(itr.ht, threadid, itr.nthreads,
                    itr.func, threadid, itr.arg);
}

static inline void hash_table_iterate(const HashTable *ht, size_t nthreads,
                                      bool (*func)(hkey_t _h, size_t threadid,
                                                   void *_arg),
                                      void *arg)
{
  ctx_assert(nthreads > 0);
  HashTableIterator ht_iter = {.ht = ht, .nthreads = nthreads,
                               .func = func, .arg = arg};

  util_multi_thread(&ht_iter, nthreads, _hash_table_iterate);
}

#endif /* HASH_TABLE_H_ */
