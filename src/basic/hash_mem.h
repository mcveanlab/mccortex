#ifndef HASH_MEM_H_
#define HASH_MEM_H_

#define REHASH_LIMIT 20
#define IDEAL_OCCUPANCY 0.75f
#define WARN_OCCUPANCY 0.9f
// bucket size must be <256
#define MAX_BUCKET_SIZE 48

// Hash table capacity is x*(2^y) where x and y are parameters
// memory is x*(2^y)*sizeof(BinaryKmer) + (2^y) * 2
static inline size_t ht_mem(size_t bktsize, size_t nbkts, size_t nbits) {
  return (bktsize * nbkts * nbits)/8 + (nbkts) * sizeof(uint8_t[2]);
}

// Returns capacity of a hash table that holds at least nkmers
size_t hash_table_cap(uint64_t nkmers, uint64_t *num_bkts_ptr, uint8_t *bkt_size_ptr);

// Returns memory required to hold nkmers
size_t hash_table_mem(uint64_t nkmers, size_t entrybits, uint64_t *nkmers_ptr);

// Returns memory used for hashtable no more than some memory limit
size_t hash_table_mem_limit(size_t memlimit, size_t entrybits, uint64_t *nkmers_ptr);

#endif /* HASH_MEM_H_ */
