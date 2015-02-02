#include "global.h"
#include "hash_table.h"
#include "hash_mem.h"
#include "util.h"

// bit macros from BitArray library used for spinlocking
#include "bit_array/bit_macros.h"

// Hash table prefetching doesn't appear to be faster
#define HASH_PREFETCH 1

static const BinaryKmer unset_bkmer = {.b = {UNSET_BKMER_WORD}};

#define ht_bckt_ptr(ht,bckt) ((ht)->table + (size_t)bckt * (ht)->bucket_size)

void hash_table_alloc(HashTable *ht, uint64_t req_capacity)
{
  uint64_t num_of_buckets, capacity;
  uint8_t bucket_size;

  capacity = hash_table_cap(req_capacity, &num_of_buckets, &bucket_size);
  uint_fast32_t hash_mask = (uint_fast32_t)(num_of_buckets - 1);

  size_t mem = capacity * sizeof(BinaryKmer) +
               num_of_buckets * sizeof(uint8_t[2]);

  char num_bkts_str[100], bkt_size_str[100], cap_str[100], mem_str[100];
  ulong_to_str(num_of_buckets, num_bkts_str);
  ulong_to_str(bucket_size, bkt_size_str);
  ulong_to_str(capacity, cap_str);
  bytes_to_str(mem, 1, mem_str);
  status("[hasht] Allocating table with %s entries, using %s", cap_str, mem_str);
  status("[hasht]  number of buckets: %s, bucket size: %s", num_bkts_str, bkt_size_str);

  // calloc is required for bucket_data to set the first element of each bucket
  // to the 0th pos
  BinaryKmer *table = ctx_malloc(capacity * sizeof(BinaryKmer));
  uint8_t (*const buckets)[2] = ctx_calloc(num_of_buckets, sizeof(uint8_t[2]));

  size_t i;
  for(i = 0; i < capacity; i++) table[i] = unset_bkmer;

  HashTable data = {
    .table = table,
    .num_of_buckets = num_of_buckets,
    .hash_mask = hash_mask,
    .bucket_size = bucket_size,
    .capacity = capacity,
    .buckets = buckets,
    .num_kmers = 0,
    .collisions = {0},
    .seed = rand()};

  memcpy(ht, &data, sizeof(data));
}

void hash_table_dealloc(HashTable *hash_table)
{
  ctx_free(hash_table->table);
  ctx_free(hash_table->buckets);
}

void hash_table_empty(HashTable *const ht)
{
  size_t i;
  BinaryKmer *table = ht->table;
  for(i = 0; i < ht->capacity; i++) table[i] = unset_bkmer;
  memset(ht->buckets, 0, ht->num_of_buckets * sizeof(uint8_t[2]));

  HashTable data = {
    .table = ht->table,
    .num_of_buckets = ht->num_of_buckets,
    .hash_mask = ht->hash_mask,
    .bucket_size = ht->bucket_size,
    .capacity = ht->capacity,
    .buckets = ht->buckets,
    .num_kmers = 0,
    .collisions = {0}};

  memcpy(ht, &data, sizeof(data));
}

static inline const BinaryKmer* hash_table_find_in_bucket_mt(const HashTable *const ht,
                                                             uint_fast32_t bucket,
                                                             const BinaryKmer bkmer)
{
  const BinaryKmer *ptr = ht_bckt_ptr(ht, bucket);
  const BinaryKmer *end = ptr + *(volatile __typeof(ht->buckets[0][0])*)&ht->buckets[bucket][HT_BSIZE];

  while(ptr < end) {
    BinaryKmer tgt = *(volatile const BinaryKmer*)ptr;
    if(binary_kmers_are_equal(bkmer, tgt)) return ptr;
    ptr++;
  }
  return NULL; // Not found
}

/*
static inline
const BinaryKmer* hash_table_find_insert_in_bucket(const HashTable *const ht,
                                                   uint_fast32_t bucket,
                                                   const BinaryKmer bkmer,
                                                   bool *found)
{
  const BinaryKmer *ptr = ht_bckt_ptr(ht, bucket);
  const BinaryKmer *end = ptr + ht->buckets[bucket][HT_BSIZE];
  const BinaryKmer *empty = NULL;

  for(; ptr < end && !binary_kmers_are_equal(bkmer, *ptr); ptr++) {
    if(!HASH_ENTRY_ASSIGNED(*ptr)) empty = ptr;
  }

  *found = (ptr < end);

  if(ptr == end && empty == NULL &&
     ht->buckets[bucket][HT_BSIZE] < ht->bucket_size)
  {
    *empty = bkmer;
    ht->num_kmers++;
    ht->buckets[bucket][HT_BITEMS]++;
    ht->buckets[bucket][HT_BSIZE]++;
  }

  return ptr < end ? ptr : empty;
}

// Code to find/insert:
// h = binary_kmer_hash(key,ht->seed+i) & ht->hash_mask;
// ptr = hash_table_find_insert_in_bucket(ht, h, key, &f);
// if(ptr != NULL) {
//   *found = f;
//   return ptr;
// }
*/

// Remember to increment ht->num_kmers
static inline BinaryKmer* hash_table_insert_in_bucket(HashTable *ht,
                                                      uint_fast32_t bucket,
                                                      const BinaryKmer bkmer)
{
  ctx_assert(ht->buckets[bucket][HT_BITEMS] < ht->bucket_size);
  BinaryKmer *ptr = ht_bckt_ptr(ht, bucket);

  if(ht->buckets[bucket][HT_BSIZE] == ht->buckets[bucket][HT_BITEMS]) {
    ptr += ht->buckets[bucket][HT_BSIZE];
    ht->buckets[bucket][HT_BSIZE]++;
  }
  else {
    // Find an entry that has been deleted from this bucket previously
    while(HASH_ENTRY_ASSIGNED(*ptr)) ptr++;
  }

  *ptr = bkmer;
  ht->buckets[bucket][HT_BITEMS]++;
  return ptr;
}

#define rehash_error_exit(ht) do { \
  ctx_msg_out = stderr; \
  hash_table_print_stats(ht); \
  die("Hash table is full"); \
} while(0)

hkey_t hash_table_find(const HashTable *const ht, const BinaryKmer key)
{
  const BinaryKmer *ptr;
  size_t i;
  uint_fast32_t h;

  #ifdef HASH_PREFETCH
    uint_fast32_t h2 = binary_kmer_hash(key,ht->seed+0) & ht->hash_mask;
    __builtin_prefetch(ht_bckt_ptr(ht, h2), 0, 1);
  #endif

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    #ifdef HASH_PREFETCH
      h = h2;
      if(ht->buckets[h][HT_BSIZE] == ht->bucket_size) {
        h2 = binary_kmer_hash(key,ht->seed+i+1) & ht->hash_mask;
        __builtin_prefetch(ht_bckt_ptr(ht, h2), 0, 1);
      }
    #else
      h = binary_kmer_hash(key,ht->seed+i) & ht->hash_mask;
    #endif

    ptr = hash_table_find_in_bucket_mt(ht, h, key);
    if(ptr != NULL) return (hkey_t)(ptr - ht->table);
    if(ht->buckets[h][HT_BSIZE] < ht->bucket_size) return HASH_NOT_FOUND;
  }

  rehash_error_exit(ht);
}

// This methods inserts an element in the next available bucket
// It doesn't check whether another element with the same key is present in the
// table used for fast loading when it is known that all the elements in the
// input have different key
hkey_t hash_table_insert(HashTable *const ht, const BinaryKmer key)
{
  const BinaryKmer *ptr;
  size_t i;
  uint_fast32_t h;
  // prefetch doesn't make sense when not searching..

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    h = binary_kmer_hash(key,ht->seed+i) & ht->hash_mask;
    if(ht->buckets[h][HT_BITEMS] < ht->bucket_size) {
      ptr = hash_table_insert_in_bucket(ht, h, key);
      ht->collisions[i]++; // only increment collisions when inserting
      ht->num_kmers++;
      return (hkey_t)(ptr - ht->table);
    }
  }

  rehash_error_exit(ht);
}

hkey_t hash_table_find_or_insert(HashTable *ht, const BinaryKmer key,
                                 bool *found)
{
  const BinaryKmer *ptr;
  size_t i;
  uint_fast32_t h;

  #ifdef HASH_PREFETCH
    uint_fast32_t h2 = binary_kmer_hash(key,ht->seed+0) & ht->hash_mask;
    __builtin_prefetch(ht_bckt_ptr(ht, h2), 0, 1);
  #endif

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    #ifdef HASH_PREFETCH
      h = h2;
      if(ht->buckets[h][HT_BSIZE] == ht->bucket_size) {
        h2 = binary_kmer_hash(key,ht->seed+i+1) & ht->hash_mask;
        __builtin_prefetch(ht_bckt_ptr(ht, h2), 0, 1);
      }
    #else
      h = binary_kmer_hash(key,ht->seed+i) & ht->hash_mask;
    #endif

    ptr = hash_table_find_in_bucket_mt(ht, h, key);

    if(ptr != NULL)  {
      *found = true;
      return (hkey_t)(ptr - ht->table);
    }
    else if(ht->buckets[h][HT_BITEMS] < ht->bucket_size) {
      *found = false;
      ptr = hash_table_insert_in_bucket(ht, h, key);
      ht->collisions[i]++; // only increment collisions when inserting
      ht->num_kmers++;
      return (hkey_t)(ptr - ht->table);
    }
  }

  rehash_error_exit(ht);
}

hkey_t hash_table_find_or_insert_mt(HashTable *ht, const BinaryKmer key,
                                    bool *found, volatile uint8_t *bktlocks)
{
  const BinaryKmer *ptr;
  size_t i;
  uint_fast32_t h;

  #ifdef HASH_PREFETCH
    uint_fast32_t h2 = binary_kmer_hash(key,ht->seed+0) & ht->hash_mask;
    __builtin_prefetch(ht_bckt_ptr(ht, h2), 0, 1);
  #endif

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    #ifdef HASH_PREFETCH
      h = h2;
      if(ht->buckets[h][HT_BSIZE] == ht->bucket_size) {
        h2 = binary_kmer_hash(key,ht->seed+i+1) & ht->hash_mask;
        __builtin_prefetch(ht_bckt_ptr(ht, h2), 0, 1);
      }
    #else
      h = binary_kmer_hash(key,ht->seed+i) & ht->hash_mask;
    #endif

    bitlock_yield_acquire(bktlocks, h);

    // We have the bucket lock so noone else can find or insert elements
    // therefore we can use non-threadsafe bucket functions
    // bitlock_acquire/release provide memory barriers
    ptr = hash_table_find_in_bucket_mt(ht, h, key);

    if(ptr != NULL)  {
      *found = true;
      bitlock_release(bktlocks, h);
      return (hkey_t)(ptr - ht->table);
    }
    else if(ht->buckets[h][HT_BITEMS] < ht->bucket_size) {
      *found = false;
      ptr = hash_table_insert_in_bucket(ht, h, key);
      bitlock_release(bktlocks, h);
      __sync_add_and_fetch((volatile uint64_t*)&ht->collisions[i], 1);
      __sync_add_and_fetch((volatile uint64_t*)&ht->num_kmers, 1);
      return (hkey_t)(ptr - ht->table);
    }
    else {
      bitlock_release(bktlocks, h);
    }
  }

  rehash_error_exit(ht);
}

// Safe to call on different entries at the same time
// NOT safe to do find() whilst doing delete()
void hash_table_delete(HashTable *const ht, hkey_t pos)
{
  uint64_t bucket = pos / ht->bucket_size;

  ctx_assert(pos != HASH_NOT_FOUND);
  ctx_assert(ht->buckets[bucket][HT_BITEMS] > 0);
  ctx_assert(ht->num_kmers > 0);
  ctx_assert(HASH_ENTRY_ASSIGNED(ht->table[pos]));

  ht->table[pos] = unset_bkmer;
  __sync_fetch_and_sub((volatile uint8_t *)&ht->buckets[bucket][HT_BITEMS], 1);
  __sync_fetch_and_sub((volatile uint64_t *)&ht->num_kmers, 1);

  ctx_assert(!HASH_ENTRY_ASSIGNED(ht->table[pos]));
}

void hash_table_print_stats_brief(const HashTable *const ht)
{
  size_t nbytes, nkeybits;
  double occupancy = (100.0 * ht->num_kmers) / ht->capacity;
  nbytes = ht->capacity * sizeof(BinaryKmer) +
           ht->num_of_buckets * sizeof(uint8_t[2]);
  nkeybits = (size_t)__builtin_ctzl(ht->num_of_buckets);

  char mem_str[50], num_buckets_str[100], num_entries_str[100], capacity_str[100];
  ulong_to_str(ht->num_of_buckets, num_buckets_str);
  bytes_to_str(nbytes, 1, mem_str);
  ulong_to_str(ht->capacity, capacity_str);
  ulong_to_str(ht->num_kmers, num_entries_str);

  status("[hash] buckets: %s [2^%zu]; bucket size: %zu; "
         "memory: %s; occupancy: %s / %s (%.2f%%)\n",
         num_buckets_str, nkeybits, (size_t)ht->bucket_size, mem_str,
         num_entries_str, capacity_str, occupancy);
}

void hash_table_print_stats(const HashTable *const ht)
{
  size_t i;
  hash_table_print_stats_brief(ht);

  if(ht->num_kmers > 0) {
    for(i = 0; i < REHASH_LIMIT; i++) {
      if(ht->collisions[i] != 0) {
        status("  collisions %2zu: %zu\n", i, (size_t)ht->collisions[i]);
      }
    }
  }
}


static inline void increment_count(hkey_t hkey, uint64_t *count)
{
  (void)hkey;
  (*count)++;
}

uint64_t hash_table_count_kmers(const HashTable *const ht)
{
  uint64_t count = 0;
  HASH_ITERATE(ht, increment_count, &count);
  return count;
}
