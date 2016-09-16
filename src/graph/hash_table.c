#include "global.h"
#include "hash_table.h"
#include "hash_mem.h"
#include "util.h"

// bit macros from BitArray library used for spinlocking
#include "bit_array/bit_macros.h"

// Hash table prefetching doesn't appear to be faster
// #define HASH_PREFETCH 1

#define ht_bckt_ptr(ht,bckt) ((ht)->table + (size_t)bckt * (ht)->bucket_size)
#define hash_table_bsize_mt(ht,bkt) (*(volatile uint8_t*)&ht->buckets[bkt][HT_BSIZE])
#define hash_table_bitems_mt(ht,bkt) (*(volatile uint8_t*)&ht->buckets[bkt][HT_BITEMS])

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
  BinaryKmer *table = ctx_calloc(capacity, sizeof(BinaryKmer));
  uint8_t (*const buckets)[2] = ctx_calloc(num_of_buckets, sizeof(uint8_t[2]));

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
  memset(ht->table, 0, ht->capacity * sizeof(BinaryKmer));
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

static inline const BinaryKmer* hash_table_find_in_bucket(const HashTable *const ht,
                                                          uint_fast32_t bucket,
                                                          BinaryKmer bkmer)
{
  const BinaryKmer *ptr = ht_bckt_ptr(ht, bucket);
  const BinaryKmer *end = ptr + hash_table_bsize(ht, bucket);
  bkmer.b[0] |= BKMER_SET_FLAG; // mark as assigned in the hash table

  while(ptr < end) {
    if(binary_kmer_eq(bkmer, *ptr)) return ptr;
    ptr++;
  }
  return NULL; // Not found
}

// Remember to increment ht->num_kmers
static inline BinaryKmer* hash_table_insert_in_bucket(HashTable *ht,
                                                      uint_fast32_t bucket,
                                                      BinaryKmer bkmer)
{
  size_t bsize = hash_table_bsize(ht, bucket);
  size_t bitems = hash_table_bitems(ht, bucket);
  ctx_assert(bitems < ht->bucket_size);
  ctx_assert(bitems <= bsize);
  BinaryKmer *ptr = ht_bckt_ptr(ht, bucket);
  bkmer.b[0] |= BKMER_SET_FLAG; // mark as assigned in the hash table

  if(bitems == bsize) {
    ptr += bsize;
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

    ptr = hash_table_find_in_bucket(ht, h, key);
    if(ptr != NULL) return (hkey_t)(ptr - ht->table);
    if(ht->buckets[h][HT_BSIZE] < ht->bucket_size) break;
  }

  return HASH_NOT_FOUND;
}

hkey_t hash_table_find_mt(HashTable *ht, const BinaryKmer key,
                          volatile uint8_t *bktlocks)
{
  const BinaryKmer *ptr;
  size_t i, bsize;
  uint_fast32_t h;

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    h = binary_kmer_hash(key,ht->seed+i) & ht->hash_mask;
    bitlock_yield_acquire(bktlocks, h);
    ptr = hash_table_find_in_bucket(ht, h, key);

    if(ptr != NULL) {
      bitlock_release(bktlocks, h);
      return (hkey_t)(ptr - ht->table);
    }

    bsize = hash_table_bsize(ht, h);
    bitlock_release(bktlocks, h);

    if(bsize < ht->bucket_size) break;
  }

  return HASH_NOT_FOUND;
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

    ptr = hash_table_find_in_bucket(ht, h, key);

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

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    h = binary_kmer_hash(key,ht->seed+i) & ht->hash_mask;
    bitlock_yield_acquire(bktlocks, h);
    ptr = hash_table_find_in_bucket(ht, h, key);

    if(ptr != NULL)  {
      *found = true;
      bitlock_release(bktlocks, h);
      return (hkey_t)(ptr - ht->table);
    }
    else if(hash_table_bitems(ht, h) < ht->bucket_size) {
      *found = false;
      ptr = hash_table_insert_in_bucket(ht, h, key);
      __sync_add_and_fetch((volatile uint64_t*)&ht->collisions[i], 1);
      __sync_add_and_fetch((volatile uint64_t*)&ht->num_kmers, 1);
      bitlock_release(bktlocks, h);
      return (hkey_t)(ptr - ht->table);
    }

    bitlock_release(bktlocks, h);
  }

  rehash_error_exit(ht);
}

// Safe to call on different entries at the same time
// NOT safe to do find() whilst doing delete()
void hash_table_delete(HashTable *const ht, hkey_t pos)
{
  uint64_t bucket = pos / ht->bucket_size, n, m;

  ctx_assert(pos != HASH_NOT_FOUND);
  ctx_assert(HASH_ENTRY_ASSIGNED(ht->table[pos]));

  memset(ht->table+pos, 0, sizeof(BinaryKmer));
  n = __sync_fetch_and_sub((volatile uint64_t *)&ht->num_kmers, 1);
  m = __sync_fetch_and_sub((volatile uint8_t *)&ht->buckets[bucket][HT_BITEMS], 1);

  ctx_assert2(n > 0, "Deleted from empty table");
  ctx_assert2(m > 0, "Deleted from empty bucket");
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

  status("[hasht] buckets: %s [2^%zu]; bucket size: %zu; ",
         num_buckets_str, nkeybits, (size_t)ht->bucket_size);
  status("[hasht] memory: %s; filled: %s / %s (%.2f%%)\n",
         mem_str, num_entries_str, capacity_str, occupancy);
}

void hash_table_print_stats(const HashTable *const ht)
{
  size_t i;
  hash_table_print_stats_brief(ht);

  if(ht->num_kmers > 0) {
    for(i = 0; i < REHASH_LIMIT; i++) {
      if(ht->collisions[i] != 0) {
        status("[hasht]  collisions %2zu: %zu\n", i, (size_t)ht->collisions[i]);
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

//
// Get hash indices sorted by kmer value
//

typedef union { BinaryKmer *bptr; hkey_t h; } BkmerPtrHkeyUnion;
static inline void _fetch_kmer_union(hkey_t hkey, const HashTable *htable,
                                     BkmerPtrHkeyUnion **bkptr)
{
  (**bkptr).bptr = htable->table + hkey;
  (*bkptr)++;
}

// Returns sorted array of hkey_t from the hash table
hkey_t* hash_table_sorted(const HashTable *htable)
{
  ctx_assert(sizeof(hkey_t) == sizeof(BinaryKmer*));
  ctx_assert(sizeof(hkey_t) == sizeof(BkmerPtrHkeyUnion));
  BkmerPtrHkeyUnion *kmers, *nxt, *end;
  nxt = kmers = ctx_malloc(sizeof(BkmerPtrHkeyUnion) * htable->num_kmers);
  end = kmers + htable->num_kmers;
  HASH_ITERATE(htable, _fetch_kmer_union, htable, &nxt);
  // Can sort ignoring that the top flag bit is set on all kmers
  qsort(kmers, htable->num_kmers, sizeof(BinaryKmer*), binary_kmers_qcmp_ptrs);
  for(nxt = kmers; nxt < end; nxt++) nxt->h = nxt->bptr - htable->table;
  return (hkey_t*)kmers;
}
