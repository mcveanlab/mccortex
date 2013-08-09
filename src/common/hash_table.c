#include "global.h"
#include "hash_table.h"
#include "util.h"

#define BSIZE 0
#define BITEMS 1

// Returns capacity
size_t hash_table_cap(size_t max_capacity_kmers,
                      uint64_t *num_bckts_ptr, uint8_t *bckt_size_ptr)
{
  uint64_t num_of_bits = 10;
  while(max_capacity_kmers / (1UL << num_of_bits) > MAX_BUCKET_SIZE) num_of_bits++;
  uint64_t num_of_buckets = 1UL << num_of_bits;
  uint64_t bucket_size = MAX2(max_capacity_kmers / num_of_buckets, 1);
  if(num_bckts_ptr != NULL) *num_bckts_ptr = num_of_buckets;
  if(bckt_size_ptr != NULL) *bckt_size_ptr = bucket_size;
  return num_of_buckets * bucket_size;
}

size_t hash_table_mem(size_t req_capacity_kmers, size_t *act_capacity_kmers)
{
  uint64_t num_of_buckets;
  size_t capacity = hash_table_cap(req_capacity_kmers, &num_of_buckets, NULL);
  if(act_capacity_kmers != NULL) *act_capacity_kmers = capacity;
  return capacity * sizeof(BinaryKmer) + num_of_buckets*sizeof(uint8_t[2]);
}

HashTable* hash_table_alloc(HashTable *htable, uint64_t req_capacity)
{
  uint64_t num_of_buckets;
  uint8_t bucket_size;
  uint64_t capacity = hash_table_cap(req_capacity, &num_of_buckets, &bucket_size);
  uint32_t hash_mask = num_of_buckets - 1;

  char capacity_str[100];
  ulong_to_str(capacity, capacity_str);
  message("[hash] Attempting to alloc table with %s entries\n", capacity_str);

  // calloc is required for bucket_data to set the first element of each bucket
  // to the 0th pos
  BinaryKmer *table = malloc2(capacity * sizeof(BinaryKmer));
  uint8_t (*const buckets)[2] = calloc2(num_of_buckets, sizeof(uint8_t[2]));

  uint64_t i, num = NUM_BKMER_WORDS * capacity, *ptr = (uint64_t*)table;
  for(i = 0; i < num; i++) ptr[i] = UNSET_BKMER;

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

static inline const BinaryKmer* hash_table_find_in_bucket(const HashTable *const htable,
                                                          uint32_t bucket,
                                                          BinaryKmer bkmer)
{
  const BinaryKmer *ptr = htable->table + (size_t)bucket * htable->bucket_size;
  const BinaryKmer *end = ptr + htable->buckets[bucket][BSIZE];

  while(ptr < end) {
    if(binary_kmers_are_equal(bkmer, *ptr)) return ptr;
    ptr++;
  }
  return NULL; // Not found
}

static inline BinaryKmer* hash_table_insert_in_bucket(HashTable *htable,
                                                      uint32_t hash,
                                                      uint32_t rehash,
                                                      BinaryKmer bkmer)
{
  assert(htable->buckets[hash][BITEMS] < htable->bucket_size);
  BinaryKmer *ptr = htable->table + (size_t)hash * htable->bucket_size;

  if(htable->buckets[hash][BSIZE] == htable->buckets[hash][BITEMS]) {
    ptr += htable->buckets[hash][BSIZE];
    htable->buckets[hash][BSIZE]++;
  }
  else {
    // Find an entry that has been deleted from this bucket previously
    while(HASH_ENTRY_ASSIGNED(*ptr)) ptr++;
  }

  *ptr = bkmer;
  htable->unique_kmers++;
  htable->buckets[hash][BITEMS]++;
  htable->collisions[rehash]++; // only increment collisions when inserting
  return ptr;
}

static inline void rehash_error_exit(const HashTable *const htable)
__attribute__((noreturn));

static inline void rehash_error_exit(const HashTable *const  htable)
{
  hash_table_print_stats(htable);
  die("Hash table is full");
}

hkey_t hash_table_find(const HashTable *const htable, const BinaryKmer const bkmer)
{
  const BinaryKmer *ptr;
  uint32_t i, hash;

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    hash = binary_kmer_hash(bkmer, i, htable->hash_mask);
    ptr = hash_table_find_in_bucket(htable, hash, bkmer);
    if(ptr != NULL) return (ptr - htable->table);
    if(htable->buckets[hash][BSIZE] < htable->bucket_size) return HASH_NOT_FOUND;
  }

  rehash_error_exit(htable);
}

// This methods inserts an element in the next available bucket
// It doesn't check whether another element with the same key is present in the
// table used for fast loading when it is known that all the elements in the
// input have different key
hkey_t hash_table_insert(HashTable *const htable, const BinaryKmer const key)
{
  uint32_t rehash, hash;
  const BinaryKmer *ptr;

  for(rehash = 0; rehash < REHASH_LIMIT; rehash++)
  {
    hash = binary_kmer_hash(key, rehash, htable->hash_mask);
    if(htable->buckets[hash][BITEMS] < htable->bucket_size) {
      ptr = hash_table_insert_in_bucket(htable, hash, rehash, key);
      return (ptr - htable->table);
    }
  }

  rehash_error_exit(htable);
}

hkey_t hash_table_find_or_insert(HashTable *htable, const BinaryKmer const key,
                                 boolean *found)
{
  const BinaryKmer *ptr;
  uint32_t rehash, hash;

  for(rehash = 0; rehash < REHASH_LIMIT; rehash++)
  {
    hash = binary_kmer_hash(key, rehash, htable->hash_mask);
    ptr = hash_table_find_in_bucket(htable, hash, key);

    if(ptr != NULL)  {
      *found = true;
      return (ptr - htable->table);
    }
    else if(htable->buckets[hash][BITEMS] < htable->bucket_size) {
      *found = false;
      ptr = hash_table_insert_in_bucket(htable, hash, rehash, key);
      return (ptr - htable->table);
    }
  }

  rehash_error_exit(htable);
}

void hash_table_delete(HashTable *const htable, hkey_t pos)
{
  BinaryKmer unset = {.b = {UNSET_BKMER}};
  htable->table[pos] = unset;
  uint64_t bucket = pos / htable->bucket_size;
  assert(htable->buckets[bucket][BITEMS] > 0);
  htable->buckets[bucket][BITEMS]--;
  htable->unique_kmers--;
}

void hash_table_print_stats(const HashTable *const htable)
{
  double occupancy = 100 * (double)htable->unique_kmers / htable->capacity;
  size_t bytes = htable->capacity * sizeof(BinaryKmer) +
                 htable->num_of_buckets * sizeof(uint8_t[2]);
  // size_t mem_height = __builtin_ctzl(htable->num_of_buckets);
  // size_t mem_width = htable->bucket_size;
  char mem_str[50], num_buckets_str[100], num_entries_str[100], capacity_str[100];
  ulong_to_str(htable->num_of_buckets, num_buckets_str);
  bytes_to_str(bytes, 1, mem_str);
  ulong_to_str(htable->capacity, capacity_str);
  ulong_to_str(htable->unique_kmers, num_entries_str);

  message("[hash table]  buckets: %s; bucket size: %zu; "
          "memory: %s; occupancy: %s / %s (%.2f%%)\n",
          num_buckets_str, (size_t)htable->bucket_size, mem_str,
          num_entries_str, capacity_str, occupancy);

  if(htable->unique_kmers > 0)
  {
    int i;
    message("  Collisions:\n");
    for(i = 0; i < REHASH_LIMIT; i++) {
      if(htable->collisions[i] != 0) {
        message("   tries %i: %zd\n", i, (size_t)htable->collisions[i]);
      }
    }
    message("\n");
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
