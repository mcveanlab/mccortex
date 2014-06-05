#include "global.h"
#include "util.h"
#include "db_graph.h"
#include "path_store.h"
#include "dna.h"
#include "binary_seq.h"

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}..
// prev = PATH_NULL if not set

void path_store_alloc(PathStore *ps, size_t mem, bool use_path_hash,
                      size_t kmers_in_hash, size_t ncols)
{
  size_t colset_bytes = roundup_bits2bytes(ncols);
  size_t pstore_mem = 0, hash_mem = 0;

  // Split memory evenly between hash and path store
  if(use_path_hash) {
    pstore_mem = mem / 2;
    hash_mem = mem / 2;
  } else {
    pstore_mem = mem;
    hash_mem = 0;
  }

  // all one block: <mem><padding><tmp><padding>
  uint8_t *block = ctx_malloc(pstore_mem + PSTORE_PADDING);

  // Paths
  PathIndex *kmer_paths = ctx_malloc(kmers_in_hash * sizeof(PathIndex));
  memset(kmer_paths, 0xff, kmers_in_hash * sizeof(PathIndex));

  char main_mem_str[100];
  bytes_to_str(pstore_mem, 1, main_mem_str);
  status("[paths] Setting up path store to use %s main", main_mem_str);

  PathStore new_paths = {.store = block,
                         .end = block + pstore_mem, .next = block,
                         .num_of_cols = ncols,
                         .colset_bytes = colset_bytes,
                         .extra_bytes = 0, .num_of_bytes = 0,
                         .num_of_paths = 0, .num_kmers_with_paths = 0,
                         .num_col_paths = 0,
                         .tmpstore = NULL, .tmpsize = 0,
                         .kmer_paths_read = kmer_paths,
                         .kmer_paths_write = kmer_paths,
                         .kmer_locks = NULL,
                         .phash = PATH_HASH_EMPTY};

  if(use_path_hash)
    path_hash_alloc(&new_paths.phash, hash_mem);

  memcpy(ps, &new_paths, sizeof(PathStore));
}

// Set up temporary memory for merging PathStores
void path_store_setup_tmp(PathStore *ps, size_t tmp_mem)
{
  size_t max_tmp_mem = (ps->end - ps->store - PSTORE_PADDING) / 2;
  ctx_assert(tmp_mem <= max_tmp_mem);

  ps->tmpsize = tmp_mem;
  ps->end -= PSTORE_PADDING + ps->tmpsize;
  ps->tmpstore = ps->end + PSTORE_PADDING;

  char mem0_str[100], mem1_str[100];
  bytes_to_str(ps->end - ps->store, 1, mem0_str);
  bytes_to_str(ps->tmpsize, 1, mem1_str);
  status("[paths] Setup tmp path memory to use %s / %s", mem0_str, mem1_str);
}

// Once tmp has been used for merging, it can be reclaimed to use generally
void path_store_release_tmp(PathStore *ps)
{
  ps->end += PSTORE_PADDING + ps->tmpsize;
  ps->tmpstore = NULL;
  ps->tmpsize = 0;

  size_t data_mem = ps->end - ps->store;

  char mem_str[100];
  bytes_to_str(data_mem, 1, mem_str);
  status("[paths] Release tmp path memory to use %s main", mem_str);
}

// Release memory
void path_store_dealloc(PathStore *ps)
{
  if(ps->kmer_paths_write != ps->kmer_paths_read) ctx_free(ps->kmer_paths_write);
  ctx_free(ps->kmer_paths_read);
  ctx_free(ps->store);
  ctx_free(ps->kmer_locks);
  path_hash_dealloc(&ps->phash);
  memset(ps, 0, sizeof(PathStore));
}

void path_store_reset(PathStore *ps, size_t nkmers_in_hash)
{
  if(ps->kmer_paths_read != NULL)
    memset(ps->kmer_paths_read, 0, nkmers_in_hash * sizeof(PathIndex));
  if(ps->kmer_paths_write != NULL && ps->kmer_paths_read != ps->kmer_paths_write)
    memset(ps->kmer_paths_write, 0, nkmers_in_hash * sizeof(PathIndex));
  ps->num_of_bytes = 0;
  ps->num_of_paths = 0;
  ps->num_kmers_with_paths = 0;
  ps->num_col_paths = 0;
  ps->next = ps->store;
  if(ps->phash.table != NULL) path_hash_reset(&ps->phash);
  // done do anything to tmpstore, kmer_locks
}

// Find a path
// returns PATH_NULL if not found, otherwise index
// path_nbytes is length in bytes of bases = (num bases + 3)/4
// query is <PathLen><PackedSeq>
PathIndex path_store_find(const PathStore *ps, PathIndex last_index,
                          const uint8_t *query, size_t path_nbytes)
{
  uint8_t *packed;
  size_t offset = sizeof(PathIndex) + ps->colset_bytes;
  size_t mem = sizeof(PathLen) + path_nbytes;

  while(last_index != PATH_NULL)
  {
    packed = ps->store + last_index;
    if(memcmp(packed+offset, query, mem) == 0) return last_index;
    last_index = packedpath_get_prev(packed);
  }

  return PATH_NULL;
}

// End up with both kmer_paths and kmer_paths_update pointing to the same array
void path_store_combine_updated_paths(PathStore *pstore)
{
  if(pstore->kmer_paths_write && !pstore->kmer_paths_read) {
    pstore->kmer_paths_read = pstore->kmer_paths_write;
  }
  else if(!pstore->kmer_paths_write && pstore->kmer_paths_read) {
    pstore->kmer_paths_write = pstore->kmer_paths_read;
  }
  else if(pstore->kmer_paths_write && pstore->kmer_paths_read &&
          pstore->kmer_paths_write != pstore->kmer_paths_read)
  {
    ctx_free(pstore->kmer_paths_read);
    pstore->kmer_paths_read = pstore->kmer_paths_write;
    status("[PathStore] Path stores merged");
  }
  ctx_assert(pstore->kmer_paths_read == pstore->kmer_paths_write);
  ctx_assert(pstore->kmer_paths_read != NULL);
}

// Add to PathHash
// packed points to <len><seq>
static inline void _path_store_add_to_hash(PathStore *ps, hkey_t hkey,
                                           PathIndex pindex,
                                           const uint8_t *packed)
{
  if(ps->phash.table != NULL) {
    size_t phash_pos;
    int pret = path_hash_find_or_insert_mt(&ps->phash, hkey, packed,
                                           ps->store, ps->colset_bytes,
                                           &phash_pos);
    if(pret == 1) { // inserted
      path_hash_set_pindex(&ps->phash, phash_pos, pindex);
    }
  }
}

// Always adds!
// Only call this function if you're sure your path is unique
// colset points to <colset>, seq points to <seq>
PathIndex path_store_add_packed(PathStore *ps, hkey_t hkey, PathIndex last_index,
                                Orientation orient, PathLen plen,
                                const uint8_t *colset, const uint8_t *seq)
{
  PathLen plen_orient = packedpath_combine_lenorient(plen, orient);

  size_t pbytes = (plen+3)/4;
  size_t mem = packedpath_mem2(ps->colset_bytes, pbytes) + ps->extra_bytes;

  if(ps->next + mem >= ps->end) die("Out of memory for paths");

  uint8_t *ptr = ps->next;

  memcpy(ptr, &last_index, sizeof(PathIndex));
  ptr += sizeof(PathIndex);
  memcpy(ptr, colset, ps->colset_bytes);
  ptr += ps->colset_bytes;
  memcpy(ptr, &plen_orient, sizeof(PathLen));
  ptr += sizeof(PathLen);
  memcpy(ptr, seq, pbytes);
  ptr += pbytes;
  // Set extra byte to max count
  if(ps->extra_bytes == 1)
    *ptr = 255;

  PathIndex pindex = (PathIndex)(ps->next - ps->store);
  ps->next += mem;
  ps->num_of_paths++;
  ps->num_kmers_with_paths += (last_index == PATH_NULL);
  ps->num_of_bytes += packedpath_mem2(ps->colset_bytes, (plen+3)/4);

  // Add to PathHash (does nothing if there is no PathHash in use)
  const uint8_t *len_seq = ptr - pbytes - sizeof(PathLen);
  _path_store_add_to_hash(ps, hkey, pindex, len_seq);

  return pindex;
}

static inline void _path_count_hist(size_t idx, uint64_t *hist,
                                    const PathStore *pstore)
{
  PathIndex pindex = pstore->phash.table[idx].pindex;
  uint8_t *ptr = pstore->store+pindex;
  uint8_t count = ptr[packedpath_mem(ptr,pstore->colset_bytes)];
  hist[count]++;
}

// Get count distribution
// `hist` must be 256*sizeof(uint64_t)
void path_store_counts_histogram(PathStore *pstore, uint64_t *hist)
{
  ctx_assert(pstore->phash.table != NULL);
  PHASH_ITERATE(&pstore->phash, _path_count_hist, hist, pstore);
}

//
// Printing functions
//

void path_store_print_status(const PathStore *pstore)
{
  char paths_str[100], mem_str[100], kmers_str[100];

  ulong_to_str(pstore->num_of_paths, paths_str);
  bytes_to_str(pstore->num_of_bytes, 1, mem_str);
  ulong_to_str(pstore->num_kmers_with_paths, kmers_str);

  status("[paths] %s paths, %s path-bytes, %s kmers",
         paths_str, mem_str, kmers_str);
}

void path_store_print_path(const PathStore *paths, PathIndex pindex)
{
  PathIndex prev;
  PathLen len;
  Orientation orient;

  prev = packedpath_get_prev(paths->store + pindex);
  packedpath_get_len_orient(paths->store+pindex, paths->colset_bytes,
                            &len, &orient);

  const uint8_t *packed = paths->store + pindex;
  const uint8_t *colbitset = packedpath_get_colset(packed);
  const uint8_t *seq = packedpath_seq(packed, paths->colset_bytes);

  size_t i;
  printf("%8zu: ", (size_t)pindex);
  if(prev == PATH_NULL) printf("    NULL");
  else printf("%8zu", (size_t)prev);
  printf(" (cols:");
  for(i = 0; i < paths->num_of_cols; i++)
    if(bitset_get(colbitset, i)) printf(" %zu", i);
  printf(")[%u]: ", len);
  for(i = 0; i < len; i++)
    putc(dna_nuc_to_char(binary_seq_get(seq, i)), stdout);
  putc('\n', stdout);
}

void path_store_print_all(const PathStore *paths)
{
  PathIndex pindex = 0, store_size = (PathIndex)(paths->next - paths->store);

  while(pindex < store_size) {
    path_store_print_path(paths, pindex);
    pindex += packedpath_mem(paths->store+pindex, paths->colset_bytes);
  }
}

// packed points to <PathLen><PackedSeq>
void print_path(hkey_t hkey, const uint8_t *packed, const PathStore *pstore)
{
  size_t i;
  PathLen plen;
  Orientation orient;
  const uint8_t *seq = packed+sizeof(PathLen);

  packedpath_get_len_orient(packed-pstore->colset_bytes-sizeof(PathIndex),
                            pstore->colset_bytes, &plen, &orient);

  // Convert packed path to string
  char nucs[plen+1];
  for(i = 0; i < plen; i++) nucs[i] = dna_nuc_to_char(binary_seq_get(seq, i));
  nucs[plen] = '\0';

  // Print
  status("Path: %zu:%i len %zu %s", (size_t)hkey, orient, (size_t)plen, nucs);
}

//
// Data checks for debugging / testing
//

// Check data if exactly filled by packed paths
bool path_store_data_integrity_check(const uint8_t *data, size_t size,
                                     size_t colbytes, size_t extra_bytes)
{
  const uint8_t *ptr = data, *end = data+size;
  PathLen plen, nbytes;
  PathIndex prev;
  while(ptr < end) {
    prev = packedpath_get_prev(ptr);
    ctx_assert_ret(prev == PATH_NULL || prev < size);

    plen = packedpath_get_len(ptr, colbytes);
    nbytes = packedpath_len_nbytes(plen);
    ctx_assert_ret2(nbytes <= size && ptr + nbytes <= end,
               "nbytes: %zu size: %zu", (size_t)nbytes, (size_t)size);

    ptr += packedpath_mem2(colbytes, nbytes) + extra_bytes;
  }
  ctx_assert_ret2(ptr == end, "data: %p end: %p ptr: %p", data, end, ptr);
  return true;
}

bool path_store_integrity_check(const PathStore *pstore)
{
  ctx_assert_ret(pstore->next >= pstore->store);
  ctx_assert_ret(pstore->next <= pstore->end);
  size_t mem = pstore->next - pstore->store;
  return path_store_data_integrity_check(pstore->store, mem,
                                         pstore->colset_bytes,
                                         pstore->extra_bytes);
}
