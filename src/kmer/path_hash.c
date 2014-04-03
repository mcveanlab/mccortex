#include "global.h"
#include "path_hash.h"
#include "hash_table.h"
#include "util.h"
#include "city.h"

// Entry is [bkmer:8][pindex:5][seq:3] :bytes

struct KPEntryStruct
{
  BinaryKmer bkmer;
  PathIndex pindex:40, len:16, seq:8;
} __attribute((packed));

void path_hash_alloc(PathHash *phash, size_t mem_in_bytes)
{
  size_t i, cap_entries; uint64_t num_bkts = 0; uint8_t bkt_size = 0;

  // Decide on hash table capacity based on how much memory we can use
  cap_entries = mem_in_bytes / sizeof(KPEntry);
  hash_table_cap(cap_entries, &num_bkts, &bkt_size);
  cap_entries = num_bkts * bkt_size;

  size_t bktlocks_mem = roundup_bits2bytes(num_bkts);
  size_t mem = cap_entries*sizeof(KPEntry) + bktlocks_mem;

  char num_bkts_str[100], bkt_size_str[100], cap_str[100], mem_str[100];
  ulong_to_str(num_bkts, num_bkts_str);
  ulong_to_str(bkt_size, bkt_size_str);
  ulong_to_str(cap_entries, cap_str);
  bytes_to_str(mem, 1, mem_str);
  status("[PathHash] Allocating table with %s entries, using %s", cap_str, mem_str);
  status("[PathHash]  number of buckets: %s bucket size: %s", num_bkts_str, bkt_size_str);

  KPEntry *table = malloc2(cap_entries * sizeof(KPEntry));
  uint8_t *bktlocks = calloc2(bktlocks_mem, sizeof(uint8_t));

  ctx_assert(num_bkts * bkt_size == cap_entries);
  ctx_assert(cap_entries > 0);
  ctx_assert(sizeof(KPEntry) == sizeof(BinaryKmer) + 8);

  for(i = 0; i < cap_entries; i++)
    table[i].bkmer.b[0] = UNSET_BKMER_WORD;

  PathHash tmp = {.table = table,
                  .num_of_buckets = num_bkts,
                  .bucket_size = bkt_size,
                  .capacity = cap_entries,
                  .mask = num_bkts - 1,
                  .num_entries = 0,
                  .bktlocks = bktlocks};

  memcpy(phash, &tmp, sizeof(PathHash));
}

void path_hash_dealloc(PathHash *phash)
{
  free((void*)phash->bktlocks);
  free(phash->table);
}

static inline bool _phash_entries_match(const KPEntry *entry,
                                        BinaryKmer bkmer,
                                        PathLen plen, const uint8_t *seq,
                                        const uint8_t *pstore, size_t colbytes)
{
  const uint8_t *seq2 = packedpath_seq(pstore+entry->pindex, colbytes);

  return (binary_kmers_are_equal(entry->bkmer, bkmer) &&
          entry->len == plen &&
          entry->seq == seq[0] &&
          (plen <= 4 || memcmp(seq, seq2, (plen+3)/4) == 0));
}

// You must acquire the lock on the kmer before adding
// packed points to <PathLen><PackedSeq>
// *pos is set to the index of the entry if inserted or found
// Returns:
//   1  inserted
//   0  found
//  -1  out of memory
// Thread Safe: uses bucket level locks
int path_hash_find_or_insert_mt(PathHash *restrict phash, BinaryKmer bkmer,
                                const uint8_t *restrict packed,
                                const uint8_t *restrict pstore, size_t colbytes,
                                size_t *restrict pos)
{
  ctx_assert(phash->table != NULL);

  const uint64_t mask = phash->mask;
  const uint8_t *seq = packed+sizeof(PathLen);
  KPEntry *entry, *end;
  PathLen plen;

  memcpy(&plen, packed, sizeof(PathLen));
  plen &= PP_LENMASK;

  size_t i, path_bytes = (plen+3)/4, mem = sizeof(PathLen) + path_bytes;
  uint64_t hash = bkmer.b[0];

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    #if NUM_BKMER_WORDS > 1
      hash = CityHash64WithSeeds((const char*)packed, mem, hash, bkmer.b[1]-i);
    #else
      hash = CityHash64WithSeeds((const char*)packed, mem, hash, i);
    #endif

    hash &= mask;
    bitlock_yield_acquire(phash->bktlocks, hash);

    entry = phash->table + hash * phash->bucket_size;
    end = entry + phash->bucket_size;

    for(; entry < end; entry++)
    {
      if(!HASH_ENTRY_ASSIGNED(entry->bkmer))
      {
        *entry = (KPEntry){.bkmer = bkmer, .pindex = 0xffffffffff,
                           .len = plen, .seq = seq[0]};
        *pos = entry - phash->table;
        bitlock_release(phash->bktlocks, hash);
        return 1;
      }
      else if(_phash_entries_match(entry, bkmer, plen, seq, pstore, colbytes))
      {
        *pos = entry - phash->table;
        bitlock_release(phash->bktlocks, hash);
        return 0;
      }
    }

    bitlock_release(phash->bktlocks, hash);
  }

  return -1; // Out of space
}

void path_hash_set_pindex(PathHash *phash, size_t pos, PathIndex pindex)
{
  phash->table[pos].pindex = pindex;
}

// Get pindex of a path
PathIndex path_hash_get_pindex(PathHash *phash, size_t pos)
{
  return phash->table[pos].pindex;
}
