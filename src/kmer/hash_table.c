#include "global.h"
#include "hash_table.h"
#include "util.h"

#define BSIZE 0
#define BITEMS 1

// Hash table prefetching doesn't appear to be faster
#define HASH_PREFETCH 1

static const BinaryKmer unset_bkmer = {.b = {UNSET_BKMER_WORD}};

#define ht_bckt_ptr(ht,bckt) ((ht)->table + (size_t)bckt * (ht)->bucket_size)

// Returns capacity
size_t hash_table_cap(size_t nkmers, boolean above_nkmers,
                      uint64_t *num_bckts_ptr, uint8_t *bckt_size_ptr)
{
  assert(nkmers < HASH_NOT_FOUND);
  uint64_t num_of_bits = 10;
  while(nkmers / (1UL << num_of_bits) > MAX_BUCKET_SIZE) num_of_bits++;
  uint64_t num_of_buckets = 1UL << num_of_bits;
  if(above_nkmers) nkmers += num_of_buckets-1;
  uint64_t bucket_size = MAX2(nkmers / num_of_buckets, 1);
  if(num_bckts_ptr != NULL) *num_bckts_ptr = num_of_buckets;
  if(bckt_size_ptr != NULL) *bckt_size_ptr = (uint8_t)bucket_size;
  return num_of_buckets * bucket_size;
}

size_t hash_table_mem(size_t nkmers, boolean above_nkmers, size_t *act_capacty_kmers)
{
  uint64_t num_of_buckets;
  size_t capacity = hash_table_cap(nkmers, above_nkmers, &num_of_buckets, NULL);
  if(act_capacty_kmers != NULL) *act_capacty_kmers = capacity;
  return capacity * sizeof(BinaryKmer) + num_of_buckets*sizeof(uint8_t[2]);
}

HashTable* hash_table_alloc(HashTable *htable, uint64_t req_capacity)
{
  uint64_t num_of_buckets;
  uint8_t bucket_size;
  uint64_t capacity = hash_table_cap(req_capacity, true, &num_of_buckets, &bucket_size);
  uint_fast32_t hash_mask = (uint_fast32_t)(num_of_buckets - 1);

  char capacity_str[100];
  ulong_to_str(capacity, capacity_str);
  status("[hash] Attempting to alloc table with %s entries\n", capacity_str);

  // calloc is required for bucket_data to set the first element of each bucket
  // to the 0th pos
  BinaryKmer *table = malloc2(capacity * sizeof(BinaryKmer));
  uint8_t (*const buckets)[2] = calloc2(num_of_buckets, sizeof(uint8_t[2]));

  size_t i;
  for(i = 0; i < capacity; i++) table[i] = unset_bkmer;

  HashTable data = {
    .table = table,
    .num_of_buckets = num_of_buckets,
    .hash_mask = hash_mask,
    .bucket_size = bucket_size,
    .capacity = capacity,
    .buckets = buckets,
    .unique_kmers = 0,
    .collisions = {0}};

  memcpy(htable, &data, sizeof(data));

  return htable;
}

void hash_table_dealloc(HashTable *hash_table)
{
  free(hash_table->table);
  free(hash_table->buckets);
}

void hash_table_empty(HashTable *const htable)
{
  size_t i;
  BinaryKmer *table = htable->table;
  for(i = 0; i < htable->capacity; i++) table[i] = unset_bkmer;
  memset(htable->buckets, 0, htable->num_of_buckets * sizeof(uint8_t[2]));

  HashTable data = {
    .table = htable->table,
    .num_of_buckets = htable->num_of_buckets,
    .hash_mask = htable->hash_mask,
    .bucket_size = htable->bucket_size,
    .capacity = htable->capacity,
    .buckets = htable->buckets,
    .unique_kmers = 0,
    .collisions = {0}};

  memcpy(htable, &data, sizeof(data));
}

static inline const BinaryKmer* hash_table_find_in_bucket(const HashTable *const htable,
                                                          uint_fast32_t bucket,
                                                          const BinaryKmer bkmer)
{
  const BinaryKmer *ptr = ht_bckt_ptr(htable, bucket);
  const BinaryKmer *end = ptr + htable->buckets[bucket][BSIZE];

  while(ptr < end) {
    if(binary_kmers_are_equal(bkmer, *ptr)) return ptr;
    ptr++;
  }
  return NULL; // Not found
}

/*
static inline
const BinaryKmer* hash_table_find_insert_in_bucket(const HashTable *const htable,
                                                   uint_fast32_t bucket,
                                                   const BinaryKmer bkmer,
                                                   boolean *found)
{
  const BinaryKmer *ptr = ht_bckt_ptr(htable, bucket);
  const BinaryKmer *end = ptr + htable->buckets[bucket][BSIZE];
  const BinaryKmer *empty = NULL;

  for(; ptr < end && !binary_kmers_are_equal(bkmer, *ptr); ptr++) {
    if(!HASH_ENTRY_ASSIGNED(*ptr)) empty = ptr;
  }

  *found = (ptr < end);

  if(ptr == end && empty == NULL &&
     htable->buckets[bucket][BSIZE] < htable->bucket_size)
  {
    *empty = bkmer;
    htable->unique_kmers++;
    htable->buckets[bucket][BITEMS]++;
    htable->buckets[bucket][BSIZE]++;
  }

  return ptr < end ? ptr : empty;
}

// Code to find/insert:
// h = binary_kmer_hash(key,i) & htable->hash_mask;
// ptr = hash_table_find_insert_in_bucket(htable, h, key, &f);
// if(ptr != NULL) {
//   *found = f;
//   return ptr;
// }
*/

static inline BinaryKmer* hash_table_insert_in_bucket(HashTable *htable,
                                                      uint_fast32_t bucket,
                                                      const BinaryKmer bkmer)
{
  assert(htable->buckets[bucket][BITEMS] < htable->bucket_size);
  BinaryKmer *ptr = ht_bckt_ptr(htable, bucket);

  if(htable->buckets[bucket][BSIZE] == htable->buckets[bucket][BITEMS]) {
    ptr += htable->buckets[bucket][BSIZE];
    htable->buckets[bucket][BSIZE]++;
  }
  else {
    // Find an entry that has been deleted from this bucket previously
    while(HASH_ENTRY_ASSIGNED(*ptr)) ptr++;
  }

  *ptr = bkmer;
  htable->unique_kmers++;
  htable->buckets[bucket][BITEMS]++;
  return ptr;
}

// static inline void rehash_error_exit(const HashTable *const htable)
// __attribute__((noreturn));

// static inline void rehash_error_exit(const HashTable *const htable)
// {
//   hash_table_print_stats(htable);
//   die("Hash table is full");
// }

#define rehash_error_exit(ht) \
        { hash_table_print_stats(htable); die("Hash table is full"); }

hkey_t hash_table_find(const HashTable *const htable, const BinaryKmer key)
{
  const BinaryKmer *ptr;
  size_t i;
  uint_fast32_t h;

  #ifdef HASH_PREFETCH
    uint_fast32_t h2 = binary_kmer_hash(key,0) & htable->hash_mask;
    __builtin_prefetch(ht_bckt_ptr(htable, h2), 0, 1);
  #endif

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    #ifdef HASH_PREFETCH
      h = h2;
      if(htable->buckets[h][BSIZE] == htable->bucket_size) {
        h2 = binary_kmer_hash(key,i+1) & htable->hash_mask;
        __builtin_prefetch(ht_bckt_ptr(htable, h2), 0, 1);
      }
    #else
      h = binary_kmer_hash(key,i) & htable->hash_mask;
    #endif

    ptr = hash_table_find_in_bucket(htable, h, key);
    if(ptr != NULL) return (hkey_t)(ptr - htable->table);
    if(htable->buckets[h][BSIZE] < htable->bucket_size) return HASH_NOT_FOUND;
  }

  rehash_error_exit(htable);
}

// This methods inserts an element in the next available bucket
// It doesn't check whether another element with the same key is present in the
// table used for fast loading when it is known that all the elements in the
// input have different key
hkey_t hash_table_insert(HashTable *const htable, const BinaryKmer key)
{
  const BinaryKmer *ptr;
  size_t i;
  uint_fast32_t h;
  // prefetch doesn't make sense when not searching..

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    h = binary_kmer_hash(key,i) & htable->hash_mask;
    if(htable->buckets[h][BITEMS] < htable->bucket_size) {
      ptr = hash_table_insert_in_bucket(htable, h, key);
      htable->collisions[i]++; // only increment collisions when inserting
      return (hkey_t)(ptr - htable->table);
    }
  }

  rehash_error_exit(htable);
}

hkey_t hash_table_find_or_insert(HashTable *htable, const BinaryKmer key,
                                 boolean *found)
{
  const BinaryKmer *ptr;
  size_t i;
  uint_fast32_t h;

  #ifdef HASH_PREFETCH
    uint_fast32_t h2 = binary_kmer_hash(key,0) & htable->hash_mask;
    __builtin_prefetch(ht_bckt_ptr(htable, h2), 0, 1);
  #endif

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    #ifdef HASH_PREFETCH
      h = h2;
      if(htable->buckets[h][BSIZE] == htable->bucket_size) {
        h2 = binary_kmer_hash(key,i+1) & htable->hash_mask;
        __builtin_prefetch(ht_bckt_ptr(htable, h2), 0, 1);
      }
    #else
      h = binary_kmer_hash(key,i) & htable->hash_mask;
    #endif

    ptr = hash_table_find_in_bucket(htable, h, key);

    if(ptr != NULL)  {
      *found = true;
      return (hkey_t)(ptr - htable->table);
    }
    else if(htable->buckets[h][BITEMS] < htable->bucket_size) {
      *found = false;
      ptr = hash_table_insert_in_bucket(htable, h, key);
      htable->collisions[i]++; // only increment collisions when inserting
      return (hkey_t)(ptr - htable->table);
    }
  }

  rehash_error_exit(htable);
}

void hash_table_delete(HashTable *const htable, hkey_t pos)
{
  uint64_t bucket = pos / htable->bucket_size;

  assert(pos != HASH_NOT_FOUND);
  assert(htable->buckets[bucket][BITEMS] > 0);
  assert(htable->unique_kmers > 0);
  assert(HASH_ENTRY_ASSIGNED(htable->table[pos]));

  htable->table[pos] = unset_bkmer;
  htable->buckets[bucket][BITEMS]--;
  htable->unique_kmers--;
}

void hash_table_print_stats(const HashTable *const htable)
{
  size_t i, nbytes, nkeybits;
  double occupancy = (100.0 * htable->unique_kmers) / htable->capacity;
  nbytes = htable->capacity * sizeof(BinaryKmer) +
           htable->num_of_buckets * sizeof(uint8_t[2]);
  nkeybits = (size_t)__builtin_ctzl(htable->num_of_buckets);

  char mem_str[50], num_buckets_str[100], num_entries_str[100], capacity_str[100];
  ulong_to_str(htable->num_of_buckets, num_buckets_str);
  bytes_to_str(nbytes, 1, mem_str);
  ulong_to_str(htable->capacity, capacity_str);
  ulong_to_str(htable->unique_kmers, num_entries_str);

  status("[hash table] buckets: %s [2^%zu]; bucket size: %zu; "
         "memory: %s; occupancy: %s / %s (%.2f%%)\n",
         num_buckets_str, nkeybits, (size_t)htable->bucket_size, mem_str,
         num_entries_str, capacity_str, occupancy);

  if(htable->unique_kmers > 0)
  {
    for(i = 0; i < REHASH_LIMIT; i++) {
      if(htable->collisions[i] != 0) {
        status("  collisions %zu: %zu\n", i, (size_t)htable->collisions[i]);
      }
    }
  }
}


static inline void increment_count(hkey_t hkey, uint64_t *count)
{
  (void)hkey;
  (*count)++;
}

uint64_t hash_table_count_assigned_nodes(const HashTable *const htable)
{
  uint64_t count = 0;
  HASH_TRAVERSE(htable, increment_count, &count);
  return count;
}
