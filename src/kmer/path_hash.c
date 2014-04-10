#include "global.h"
#include "path_hash.h"
#include "util.h"
#include "city.h"

// Entry is [hkey:5][pindex:5] = 10 bytes

// We compare with REHASH_LIMIT(16)*bucket_size(<255) = 4080
// so we need 12 bits to have 2^12 = 4096 possibilities

// (1-(1/(2^12)))^4080 = 0.369 = 37% of entries would have zero collisions
// (1-(1/(2^16)))^4080 = 0.939 = 94% of entries would have zero collisions

// Packed structure is 10 bytes
// Do not use pointes to fields in this struct - they are not aligned
struct KPEntryStruct
{
  // 5 bytes each
  hkey_t hkey:40;
  PathIndex pindex:40;
} __attribute((packed));

#define PATH_HASH_UNSET (0xffffffffff)
#define PATH_HASH_ENTRY_ASSIGNED(x) ((x).hkey != PATH_HASH_UNSET)

void path_hash_alloc(PathHash *phash, size_t mem_in_bytes)
{
  size_t cap_entries; uint64_t num_bkts = 0; uint8_t bkt_size = 0;

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
  status("[PathHash]  number of buckets: %s, bucket size: %s", num_bkts_str, bkt_size_str);

  KPEntry *table = malloc2(cap_entries * sizeof(KPEntry));
  uint8_t *bktlocks = calloc2(bktlocks_mem, sizeof(uint8_t));

  ctx_assert(num_bkts * bkt_size == cap_entries);
  ctx_assert(cap_entries > 0);
  ctx_assert(sizeof(KPEntry) == 10);

  // Table all set to 1 to indicate empty
  memset(table, 0xff, cap_entries * sizeof(KPEntry));

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
  ctx_free((void*)phash->bktlocks);
  ctx_free(phash->table);
}

static inline bool _phash_entries_match(const KPEntry *entry,
                                        hkey_t hkey, PathLen plen,
                                        const uint8_t *seq,
                                        const uint8_t *pstore, size_t colbytes)
{
  const uint8_t *seq2 = packedpath_seq(pstore+entry->pindex, colbytes);
  return (hkey == entry->hkey && memcmp(seq, seq2, (plen+3)/4) == 0);
}

// You must acquire the lock on the kmer before adding
// packed points to <PathLen><PackedSeq>
// *pos is set to the index of the entry if inserted or found
// Returns:
//   1  inserted
//   0  found
//  -1  out of memory
// Thread Safe: uses bucket level locks
int path_hash_find_or_insert_mt(PathHash *restrict phash, hkey_t hkey,
                                const uint8_t *restrict packed,
                                const uint8_t *restrict pstore, size_t colbytes,
                                size_t *restrict pos)
{
  ctx_assert(phash->table != NULL);
  ctx_assert(hkey < PATH_HASH_UNSET);

  const uint64_t mask = phash->mask;
  const uint8_t *seq = packed+sizeof(PathLen);
  KPEntry *entry, *end;
  PathLen plen;

  memcpy(&plen, packed, sizeof(PathLen));
  plen &= PP_LENMASK;

  size_t i, path_bytes = (plen+3)/4, mem = sizeof(PathLen) + path_bytes;
  uint64_t hash = hkey;

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    hash = CityHash64WithSeeds((const char*)packed, mem, hash, i);

    hash &= mask;
    bitlock_yield_acquire(phash->bktlocks, hash);

    entry = phash->table + hash * phash->bucket_size;
    end = entry + phash->bucket_size;

    for(; entry < end; entry++)
    {
      if(!PATH_HASH_ENTRY_ASSIGNED(*entry))
      {
        *entry = (KPEntry){.hkey = hkey, .pindex = PATH_HASH_UNSET};
        *pos = entry - phash->table;
        bitlock_release(phash->bktlocks, hash);
        return 1;
      }
      else if(_phash_entries_match(entry, hkey, plen, seq, pstore, colbytes))
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
