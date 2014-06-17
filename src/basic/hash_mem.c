#include "global.h"
#include "hash_mem.h"

// Returns capacity of a hash table that holds at least nkmers
size_t hash_table_cap(size_t nkmers, uint64_t *num_bkts_ptr, uint8_t *bkt_size_ptr)
{
  ctx_assert(nkmers < UINT64_MAX>>1);
  uint64_t num_of_buckets, bucket_size, num_of_bits = 10;
  while(nkmers / (1UL << num_of_bits) > MAX_BUCKET_SIZE) num_of_bits++;
  num_of_buckets = 1UL << num_of_bits;
  bucket_size = MAX2((nkmers+num_of_buckets-1) / num_of_buckets, 1);
  if(num_bkts_ptr != NULL) *num_bkts_ptr = num_of_buckets;
  if(bkt_size_ptr != NULL) *bkt_size_ptr = (uint8_t)bucket_size;
  return num_of_buckets * bucket_size;
}

// Returns memory required to hold nkmers
size_t hash_table_mem(size_t nkmers, size_t entrybits, size_t *nkmers_ptr)
{
  uint64_t num_of_buckets, capacity; uint8_t bktsize;
  capacity = hash_table_cap(nkmers, &num_of_buckets, &bktsize);
  if(nkmers_ptr != NULL) *nkmers_ptr = capacity;
  return ht_mem(bktsize,num_of_buckets,entrybits);
}

// Returns memory used for hashtable no more than some memory limit
size_t hash_table_mem_limit(size_t memlimit, size_t entrybits, size_t *nkmers_ptr)
{
  size_t bktsize, num_of_bits = 10, num_of_buckets = 1UL<<num_of_bits, num_of_kmers;

  while(ht_mem(MAX_BUCKET_SIZE, num_of_buckets, entrybits) < memlimit) {
    num_of_bits++;
    num_of_buckets = 1UL << num_of_bits;
  }

  bktsize = (memlimit - num_of_buckets*sizeof(uint8_t[2])) /
            ((num_of_buckets * entrybits) /8);

  if(bktsize == 0) {
    num_of_bits--;
    num_of_buckets = 1UL << num_of_bits;
    num_of_kmers = bktsize * num_of_buckets;
    bktsize = MAX2(num_of_kmers / num_of_buckets, 1);
  }

  bktsize = MIN2(bktsize, MAX_BUCKET_SIZE);

  if(nkmers_ptr != NULL) *nkmers_ptr = num_of_buckets * bktsize;

  return ht_mem(bktsize,num_of_buckets,entrybits);
}
