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
  test_status("[hash_table] test add/delete");

  HashTable ht;
  BinaryKmer bkmer0, bkmer1, bkxor = {.b = {0}}, bkresult = {.b = {0}};
  boolean found0 = false, found1 = false;
  hkey_t key0, key1;
  size_t i, t, kmers_added = 0, kmers_deleted = 0;

  hash_table_alloc(&ht, 2048);

  for(t = 0; t < NTESTS/2; t++)
  {
    // Add two bkmers then delete one
    bkmer0 = binary_kmer_random(MAX_KMER_SIZE);
    hash_table_find_or_insert(&ht, bkmer0, &found0);
    bkmer1 = binary_kmer_random(MAX_KMER_SIZE);
    hash_table_find_or_insert(&ht, bkmer1, &found1);

    kmers_added += !found0 + !found1;

    key0 = hash_table_find(&ht, bkmer0);
    TASSERT(key0 != HASH_NOT_FOUND);
    key1 = hash_table_find(&ht, bkmer1);
    TASSERT(key1 != HASH_NOT_FOUND);

    hash_table_delete(&ht, key1);
    kmers_deleted++;

    key0 = hash_table_find(&ht, bkmer0);
    TASSERT(key0 != HASH_NOT_FOUND);
    key1 = hash_table_find(&ht, bkmer1);
    TASSERT(key1 == HASH_NOT_FOUND); // we just deleted the second bkmer

    // Keep XOR of bkmers in hash
    for(i = 0; i < NUM_BKMER_WORDS; i++) bkxor.b[i] ^= bkmer0.b[i];
  }

  TASSERT(kmers_added - kmers_deleted == ht.unique_kmers);

  // Check xor of bkmers
  size_t kcount = 0;
  HASH_ITERATE(&ht, check_bkmers, &ht, &bkresult, &kcount);

  TASSERT(kcount == ht.unique_kmers);
  for(i = 0; i < NUM_BKMER_WORDS; i++) TASSERT(bkxor.b[i] == bkresult.b[i]);

  hash_table_dealloc(&ht);
}
