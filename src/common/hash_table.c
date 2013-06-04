#include "global.h"
#include "hash_table.h"

HashTable* hash_table_alloc(HashTable *htable, uint64_t req_capacity)
{
  // bucket size must be <256
  uint64_t number_bits = 10;
  while(req_capacity / (1UL << number_bits) > 128) number_bits++;

  uint64_t number_buckets = 1UL << number_bits;
  uint64_t bucket_size = (req_capacity+number_buckets-1) / number_buckets;
  uint64_t capacity = number_buckets * bucket_size;
  uint32_t hash_mask = number_buckets - 1;

  BinaryKmer *table = malloc(capacity * sizeof(BinaryKmer));

  uint64_t i, num = NUM_BITFIELDS_IN_BKMER * capacity, *ptr = (uint64_t*)table;
  for(i = 0; i < num; i++) ptr[i] = UNSET_BKMER;

  // calloc is required to set the first element of each bucket to the 0th pos
  uint8_t *bucket_data = calloc(number_buckets*2, sizeof(uint8_t));

  uint8_t *bucket_length = bucket_data;
  uint8_t *bucket_fill = bucket_data + number_buckets*sizeof(uint8_t);

  if(table == NULL || bucket_data == NULL)
  {
    if(table != NULL) free(table);
    if(bucket_data != NULL) free(bucket_data);
    warn("could not allocate hash table with %zu entries\n", (size_t)capacity);
    return NULL;
  }

  HashTable data = {
    .table = table,
    .number_buckets = number_buckets,
    .hash_mask = hash_mask,
    .bucket_size = bucket_size,
    .capacity = capacity,
    .bucket_length = bucket_length,
    .bucket_fill = bucket_fill,
    .unique_kmers = 0,
    .collisions = {0}};

  memcpy(htable, &data, sizeof(data));

  return htable;
}

void hash_table_dealloc(HashTable *hash_table)
{
  free(hash_table->table);
  free(hash_table->bucket_length);
}

static inline BinaryKmer* hash_table_find_in_bucket(const HashTable *const htable,
                                                    const uint32_t bucket,
                                                    const BinaryKmer const bkmer)
{
  BinaryKmer *ptr = htable->table + (size_t)bucket * htable->bucket_size;
  BinaryKmer *end = ptr + htable->bucket_length[bucket];

  while(ptr < end) {
    if(binary_kmers_are_equal(bkmer, *ptr)) return ptr;
    ptr++;
  }
  return NULL; // Not found
}

static inline BinaryKmer* hash_table_insert_in_bucket(HashTable *htable,
                                                      const uint32_t hash,
                                                      const uint32_t rehash,
                                                      const BinaryKmer bkmer)
{
  assert(htable->bucket_fill[hash] < htable->bucket_size);
  BinaryKmer *ptr = htable->table + (size_t)hash * htable->bucket_size;

  if(htable->bucket_length[hash] == htable->bucket_fill[hash]) {
    ptr += htable->bucket_length[hash];
    htable->bucket_length[hash]++;
  }
  else {
    // Entries deleted from this bucket
    while(HASH_ENTRY_ASSIGNED(*ptr)) ptr++;
  }

  memcpy(ptr, bkmer, sizeof(BinaryKmer));
  htable->unique_kmers++;
  htable->bucket_fill[hash]++;
  htable->collisions[rehash]++; // only increment collisions when inserting
  return ptr;
}

static inline void rehash_error_exit()
__attribute__((noreturn));

static inline void rehash_error_exit()
{
  die("Dear user - you have not allocated enough memory to contain your "
      "sequence data. Either allocate more memory (have you done your "
      "calculations right? have you allowed for sequencing errors?), or "
      "threshold more harshly on quality score, and try again. "
      "Aborting mission.\n");
}

hkey_t hash_table_find(const HashTable *const htable, const BinaryKmer const bkmer)
{
  BinaryKmer *ptr;
  uint32_t i, hash;

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    hash = binary_kmer_hash(bkmer, i, htable->hash_mask);
    ptr = hash_table_find_in_bucket(htable, hash, bkmer);
    if(ptr != NULL) return (ptr - htable->table);
    if(htable->bucket_length[hash] < htable->bucket_size) return HASH_NOT_FOUND;
  }

  rehash_error_exit();
}

// This methods inserts an element in the next available bucket
// It doesn't check whether another element with the same key is present in the
// table used for fast loading when it is known that all the elements in the
// input have different key
hkey_t hash_table_insert(HashTable *const htable, const BinaryKmer const key)
{
  uint32_t rehash, hash;

  for(rehash = 0; rehash < REHASH_LIMIT; rehash++)
  {
    hash = binary_kmer_hash(key, rehash, htable->hash_mask);
    if(htable->bucket_fill[hash] < htable->bucket_size) {
      BinaryKmer *ptr = hash_table_insert_in_bucket(htable, hash, rehash, key);
      return (ptr - htable->table);
    }
  }

  rehash_error_exit();
}

hkey_t hash_table_find_or_insert(HashTable *htable, const BinaryKmer const key,
                                 boolean *found)
{
  BinaryKmer *ptr;
  uint32_t rehash, hash;

  for(rehash = 0; rehash < REHASH_LIMIT; rehash++)
  {
    hash = binary_kmer_hash(key, rehash, htable->hash_mask);
    ptr = hash_table_find_in_bucket(htable, hash, key);

    if(ptr != NULL)  {
      *found = true;
      return (ptr - htable->table);
    }
    else if(htable->bucket_length[hash] < htable->bucket_size) {
      *found = false;
      ptr = hash_table_insert_in_bucket(htable, hash, rehash, key);
      return (ptr - htable->table);
    }
  }

  rehash_error_exit();
}

void hash_table_delete(HashTable *const htable, hkey_t pos)
{
  assert(htable->bucket_fill[bucket] > 0);
  BinaryKmer unset = {UNSET_BKMER};
  memcpy(htable->table + pos, unset, sizeof(BinaryKmer));
  uint64_t bucket = pos / htable->bucket_size;
  htable->bucket_fill[bucket]--;
}

void hash_table_print_stats(const HashTable *const htable)
{
  message("Hash table collisions:\n");
  int i;
  for(i = 0; i < REHASH_LIMIT; i++) {
    if(htable->collisions[i] != 0) {
      message("  tries %i: %zd\n", i, (size_t)htable->collisions[i]);
    }
  }
  // Print 
  message("  %.2f%% occupancy",
          100 * (double)htable->unique_kmers / htable->capacity);
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
