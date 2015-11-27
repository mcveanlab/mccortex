#include "global.h"
#include "all_tests.h"
#include "hash_table.h"
#include "binary_kmer.h"

#define NTESTS 1024

static void xor_bkmers(hkey_t key, HashTable *ht, BinaryKmer *ptr, size_t *c)
{
  size_t i;
  BinaryKmer bkmer = ht->table[key];
  for(i = 0; i < NUM_BKMER_WORDS; i++) ptr->b[i] ^= bkmer.b[i];
  (*c)++;
}

static void test_add_remove()
{
  test_status("Test add/delete to hash_table");

  HashTable ht;
  BinaryKmer bkmer0, bkmer1, bkey0, bkey1;
  BinaryKmer bkxor = {.b = {0}}, bkresult = {.b = {0}};
  bool found0 = false, found1 = false;
  hkey_t hkey0, hkey1;
  size_t i, t, kmers_added = 0, kmers_deleted = 0;
  size_t kmer_size = MAX_KMER_SIZE;

  hash_table_alloc(&ht, 2048);

  for(t = 0; t < NTESTS/2; t++)
  {
    // Add two bkmers then delete one
    bkmer0 = binary_kmer_random(kmer_size);
    bkey0 = binary_kmer_get_key(bkmer0, kmer_size);
    hash_table_find_or_insert(&ht, bkey0, &found0);

    bkmer1 = binary_kmer_random(kmer_size);
    bkey1 = binary_kmer_get_key(bkmer1, kmer_size);
    hash_table_find_or_insert(&ht, bkey1, &found1);

    kmers_added += !found0 + !found1;

    hkey0 = hash_table_find(&ht, bkey0);
    TASSERT(hkey0 != HASH_NOT_FOUND);
    hkey1 = hash_table_find(&ht, bkey1);
    TASSERT(hkey1 != HASH_NOT_FOUND);

    hash_table_delete(&ht, hkey1);
    kmers_deleted++;

    hkey0 = hash_table_find(&ht, bkey0);
    TASSERT(hkey0 != HASH_NOT_FOUND);
    hkey1 = hash_table_find(&ht, bkey1);
    TASSERT(hkey1 == HASH_NOT_FOUND); // we just deleted the second bkmer

    // Keep XOR of bkmers in hash
    for(i = 0; i < NUM_BKMER_WORDS; i++) bkxor.b[i] ^= bkey0.b[i];
  }

  TASSERT(kmers_added - kmers_deleted == ht.num_kmers);

  // Check xor of bkmers
  size_t kcount = 0;
  HASH_ITERATE(&ht, xor_bkmers, &ht, &bkresult, &kcount);

  TASSERT(kcount == ht.num_kmers);
  TASSERT(binary_kmers_are_equal(bkxor, bkresult));

  hash_table_dealloc(&ht);
}

typedef struct {
  HashTable ht;
  uint8_t *bktlocks;
  BinaryKmer *bkmers;
  size_t *nadded;
  size_t n;
} BKmerTestSet;

void load_bset(void *arg, size_t threadid)
{
  (void)threadid;
  BKmerTestSet *bset = (BKmerTestSet*)arg;
  size_t i, start = rand() % bset->n;
  bool found = false;
  for(i = start; i < bset->n; i++) {
    hash_table_find_or_insert_mt(&bset->ht, bset->bkmers[i], &found, bset->bktlocks);
    __sync_fetch_and_add((volatile size_t*)&bset->nadded[i], !found);
  }
  sched_yield(); // release the CPU
  for(i = 0; i < start; i++) {
    hash_table_find_or_insert_mt(&bset->ht, bset->bkmers[i], &found, bset->bktlocks);
    __sync_fetch_and_add((volatile size_t*)&bset->nadded[i], !found);
  }
}

static void test_hash_table_mt()
{
  // Generate 2000 random binary kmers
  // start 20 threads adding them to the hash table
  size_t i, kmer_size = MAX_KMER_SIZE;
  size_t nthreads = (rand() % 50)+1, nkmers = 1000000;

  test_status("Testing hash table multithreading %zu threads, %zu kmers",
              nthreads, nkmers);

  BKmerTestSet bset;
  bset.n = nkmers;
  hash_table_alloc(&bset.ht, bset.n*1.5);
  bset.bkmers = ctx_calloc(bset.n, sizeof(bset.bkmers[0]));
  bset.nadded = ctx_calloc(bset.n, sizeof(bset.nadded[0]));
  bset.bktlocks = ctx_calloc((bset.ht.capacity+7)/8, 1);

  for(i = 0; i < bset.n; i++)
    bset.bkmers[i] = binary_kmer_random(kmer_size);

  util_multi_thread(&bset, nthreads, load_bset);

  // Check each kmers was added exactly once
  for(i = 0; i < bset.n; i++)
    TASSERT2(bset.nadded[i] == 1, "%zu", bset.nadded[i]);

  TASSERT(bset.ht.num_kmers == nkmers);

  ctx_free(bset.bktlocks);
  ctx_free(bset.nadded);
  ctx_free(bset.bkmers);
  hash_table_dealloc(&bset.ht);
}

void test_hash_table()
{
  test_add_remove();
  test_hash_table_mt();
}
