#include "global.h"
#include "all_tests.h"
#include "hash_table.h"
#include "binary_kmer.h"

#define NTESTS 1024

static void check_bkmers(hkey_t key, HashTable *ht, BinaryKmer *ptr, size_t *c)
{
  size_t i;
  BinaryKmer bkmer = ht->table[key];
  for(i = 0; i < NUM_BKMER_WORDS; i++) ptr->b[i] ^= bkmer.b[i];
  (*c)++;
}

void test_hash_table()
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
  HASH_ITERATE(&ht, check_bkmers, &ht, &bkresult, &kcount);

  TASSERT(kcount == ht.num_kmers);
  TASSERT(binary_kmers_are_equal(bkxor, bkresult));

  hash_table_dealloc(&ht);
}
