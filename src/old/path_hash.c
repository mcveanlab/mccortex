#include "global.h"
#include "path_hash.h"
#include "util.h"
#include "misc/city.h"

// Entry is [hkey:5][pindex:5] = 10 bytes

// We compare with REHASH_LIMIT(16)*bucket_size(<255) = 4080
// so we need 12 bits to have 2^12 = 4096 possibilities

// (1-(1/(2^12)))^4080 = 0.369 = 37% of entries would have zero collisions
// (1-(1/(2^16)))^4080 = 0.939 = 94% of entries would have zero collisions

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

  KPEntry *table = ctx_malloc(cap_entries * sizeof(KPEntry));
  uint8_t *bktlocks = ctx_calloc(bktlocks_mem, sizeof(uint8_t));
  uint8_t *bucket_nitems = ctx_calloc(num_bkts, sizeof(uint8_t));

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
                  .bucket_nitems = bucket_nitems,
                  .bktlocks = bktlocks};

  memcpy(phash, &tmp, sizeof(PathHash));
}

void path_hash_dealloc(PathHash *phash)
{
  ctx_free(phash->bucket_nitems);
  ctx_free(phash->bktlocks);
  ctx_free(phash->table);
  memset(phash, 0, sizeof(PathHash));
}

void path_hash_reset(PathHash *phash)
{
  phash->num_entries = 0;
  memset(phash->table, 0xff, phash->capacity * sizeof(KPEntry));
}

static inline bool _phash_entries_match(KPEntry entry,
                                        hkey_t hkey, PathLen plen,
                                        const uint8_t *seq,
                                        const uint8_t *pstore, size_t colbytes)
{
  ctx_assert(entry.hkey != PATH_HASH_UNSET);
  ctx_assert(entry.pindex != PATH_HASH_UNSET);
  const uint8_t *seq2 = packedpath_seq(pstore+entry.pindex, colbytes);
  return (hkey == entry.hkey && memcmp(seq, seq2, (plen+3)/4) == 0);
}

// Use a bucket lock to find or add an entry
// Returns:
//   1  inserted
//   0  found
//  -1  not found and not inserted
static inline int _find_or_add_in_bucket(PathHash *restrict phash, uint64_t hash,
                                         hkey_t hkey, PathLen plen,
                                         const uint8_t *restrict seq,
                                         const uint8_t *restrict pstore,
                                         size_t colbytes,
                                         size_t *restrict pos)
{
  bitlock_yield_acquire(phash->bktlocks, hash);

  KPEntry *entry = phash->table + hash * phash->bucket_size;
  const KPEntry *end = entry + phash->bucket_size;

  for(; entry < end; entry++)
  {
    if(!PATH_HASH_ENTRY_ASSIGNED(*entry))
    {
      *entry = (KPEntry){.hkey = hkey, .pindex = PATH_HASH_UNSET};
      *pos = entry - phash->table;
      __sync_fetch_and_add((volatile uint8_t *)&phash->bucket_nitems[hash], 1);
      bitlock_release(phash->bktlocks, hash);
      return 1;
    }
    else if(_phash_entries_match(*entry, hkey, plen, seq, pstore, colbytes))
    {
      *pos = entry - phash->table;
      bitlock_release(phash->bktlocks, hash);
      return 0;
    }
  }

  bitlock_release(phash->bktlocks, hash);

  return -1;
}

// Lock free search bucket for match.
// We can traverse a full bucket without acquiring the lock first because
// items are added but never removed from the hash. This allows us to remove
// the use of locks and improve performance.
// Returns:
//   0  found
//  -1  not found
static inline int _find_in_bucket(PathHash *restrict phash, uint64_t hash,
                                  hkey_t hkey, PathLen plen,
                                  const uint8_t *restrict seq,
                                  const uint8_t *restrict pstore,
                                  size_t colbytes,
                                  size_t *restrict pos)
{
  const KPEntry *entry = phash->table + hash * phash->bucket_size;
  const KPEntry *end = entry + phash->bucket_size;

  for(; entry < end; entry++) {
    if(_phash_entries_match(*entry, hkey, plen, seq, pstore, colbytes)) {
      *pos = entry - phash->table;
      return 0;
    }
  }

  return -1;
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
                                const uint8_t *restrict pstore,
                                size_t colbytes,
                                size_t *restrict pos)
{
  ctx_assert(phash->table != NULL);
  ctx_assert(hkey < PATH_HASH_UNSET);

  const uint64_t mask = phash->mask;
  const uint8_t *seq = packed+sizeof(PathLen);
  PathLen plen;

  memcpy(&plen, packed, sizeof(PathLen));
  plen &= PP_LENMASK;

  size_t i, path_bytes = (plen+3)/4, mem = sizeof(PathLen) + path_bytes;
  uint64_t hash = hkey;
  int ret;

  for(i = 0; i < REHASH_LIMIT; i++)
  {
    hash = CityHash64WithSeeds((const char*)packed, mem, hash, i);
    hash &= mask;

    uint8_t bucket_fill = *(volatile uint8_t *)&phash->bucket_nitems[hash];

    if(bucket_fill < phash->bucket_size) {
      ret = _find_or_add_in_bucket(phash, hash, hkey, plen, seq,
                                   pstore, colbytes, pos);
    }
    else {
      ret = _find_in_bucket(phash, hash, hkey, plen, seq,
                            pstore, colbytes, pos);
    }

    if(ret >= 0) return ret;
  }

  return -1; // Out of space
}

void path_hash_set_pindex(PathHash *phash, size_t pos, PathIndex pindex)
{
  phash->table[pos].pindex = pindex;
}

// Get pindex of a path
PathIndex path_hash_get_pindex(const PathHash *phash, size_t pos)
{
  return phash->table[pos].pindex;
}
