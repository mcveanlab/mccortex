#include "global.h"
#include "util.h"
#include "db_graph.h"
#include "path_store.h"
#include "dna.h"

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}..
// prev = PATH_NULL if not set

void path_store_alloc(PathStore *paths, size_t size, size_t tmpsize, size_t ncols)
{
  size_t colset_bytes = roundup_bits2bytes(ncols);

  // all one block: <mem><padding><tmp><padding>
  uint8_t *mem = malloc2(size + PSTORE_PADDING + tmpsize + PSTORE_PADDING);
  uint8_t *tmp = tmpsize ? mem + size + PSTORE_PADDING : NULL;

  char main_mem_str[100], tmp_mem_str[100];
  bytes_to_str(size, 1, main_mem_str);
  bytes_to_str(tmpsize, 1, tmp_mem_str);
  status("[paths] Setting up path store to use %s main / %s tmp",
         main_mem_str, tmp_mem_str);

  PathStore new_paths = {.store = mem, .end = mem + size,
                         .size = size, .next = mem,
                         .num_of_cols = ncols,
                         .colset_bytes = colset_bytes,
                         .num_of_paths = 0, .num_kmers_with_paths = 0,
                         .tmpdata = tmp, .tmpsize = tmpsize};

  memcpy(paths, &new_paths, sizeof(PathStore));
}

// Once tmp has been used for merging, it can be reclaimed to use generally
void path_store_reclaim_tmp(PathStore *paths)
{
  size_t size = paths->size + PSTORE_PADDING + paths->tmpsize;

  char memstr[100];
  bytes_to_str(size, 1, memstr);
  status("[paths] Reclaiming tmp path memory to use %s main", memstr);

  PathStore new_paths = {.store = paths->store, .end = paths->store + size,
                         .size = size, .next = paths->next,
                         .num_of_cols = paths->num_of_cols,
                         .colset_bytes = paths->colset_bytes,
                         .num_of_paths = paths->num_of_paths,
                         .num_kmers_with_paths = paths->num_kmers_with_paths,
                         .tmpdata = NULL, .tmpsize = 0};

  memcpy(paths, &new_paths, sizeof(PathStore));
}

// Release memory
void path_store_dealloc(PathStore *paths)
{
  free(paths->store);
}

// Find a path
// returns PATH_NULL if not found, otherwise index
// path_nbytes is length in bytes of bases = (num bases + 3)/4
// query is <PathLen><PackedSeq>
PathIndex path_store_find(const PathStore *paths, PathIndex last_index,
                          const uint8_t *query, size_t path_nbytes)
{
  uint8_t *packed;
  size_t offset = sizeof(PathIndex) + paths->colset_bytes;
  size_t mem = sizeof(PathLen) + path_nbytes;

  while(last_index != PATH_NULL)
  {
    packed = paths->store + last_index;
    if(memcmp(packed+offset, query, mem) == 0) return last_index;
    last_index = packedpath_get_prev(packed);
  }

  return PATH_NULL;
}

// Always adds!
// Only call this function if you're sure your path is unique
PathIndex path_store_add_packed(PathStore *store, PathIndex last_index,
                                const uint8_t *packed, size_t path_nbytes)
{
  // Not already in the PathStore
  size_t mem = packedpath_mem2(store->colset_bytes, path_nbytes);

  // Copy path (may already be in place)
  if(packed != store->next)
  {
    if(store->next + mem >= store->end) die("Out of memory for paths");

    uint8_t *ptr = store->next;
    memcpy(ptr, &last_index, sizeof(PathIndex));
    ptr += sizeof(PathIndex);
    memcpy(ptr, packed, store->colset_bytes+sizeof(PathLen)+path_nbytes);
  }

  PathIndex pindex = (PathIndex)(store->next - store->store);
  store->next += mem;
  store->num_of_paths++;
  store->num_kmers_with_paths += (last_index == PATH_NULL);
  return pindex;
}

// Find or add a path into the PathStore
// last_index is index of the last path belonging to the kmer which owns the
// new path that is to be inserted
// Returns match PathIndex if found, otherwise PATH_NULL
PathIndex path_store_find_or_add_packed(PathStore *paths, PathIndex last_index,
                                        const uint8_t *packed, size_t path_nbytes,
                                        boolean *inserted)
{
  size_t i;
  // query points to <PathLen><PackedSeq>
  const uint8_t *query = packed + sizeof(PathIndex) + paths->colset_bytes;
  PathIndex match = path_store_find(paths, last_index, query, path_nbytes);

  if(match == PATH_NULL)
  {
    match = path_store_add_packed(paths, last_index, packed, path_nbytes);
    *inserted = true;
  }
  else {
    // Already in path store, just update colour bitset
    uint8_t *dst = packedpath_get_colset(paths->store + match);
    const uint8_t *src = packedpath_get_colset(packed);
    for(i = 0; i < paths->colset_bytes; i++) dst[i] |= src[i];
    *inserted = false;
  }

  return match;
}

// Specify if we should try to find a duplicate first
static inline
PathIndex _path_store_find_or_add_packed(PathStore *store, PathIndex last_index,
                                         const uint8_t *packed, size_t path_nbytes,
                                         boolean find, boolean *added)
{
  if(find) {
    return path_store_find_or_add_packed(store, last_index, packed,
                                         path_nbytes, added);
  } else {
    *added = true;
    return path_store_add_packed(store, last_index, packed, path_nbytes);
  }
}

// Add a PackedPath, using a FileFilter to reduce to a subset of colours
// Returns PATH_NULL if no colours set in colour subset
PathIndex path_store_find_or_add_packed2(PathStore *store, PathIndex last_index,
                                         const uint8_t *packed, size_t path_nbytes,
                                         const FileFilter *fltr, boolean find,
                                         boolean *added)
{
  size_t i, intocol, fromcol, tmp;
  const size_t packed_bitset_bytes = roundup_bits2bytes(fltr->filencols);

  if(path_store_fltr_compatible(store,fltr)) {
    return _path_store_find_or_add_packed(store, last_index, packed,
                                         path_nbytes, find, added);
  }

  const uint8_t *packed_bitset = packed + sizeof(PathIndex);
  uint8_t *store_bitset = store->next + sizeof(PathIndex);

  // Check we have enough memory to add
  size_t mem = packedpath_mem2(store->colset_bytes, path_nbytes);
  if(store->next + mem >= store->end) die("Out of memory for paths");

  // Write
  memcpy(store->next, &last_index, sizeof(PathIndex));

  // Clear memory for colour bitset
  memset(store_bitset, 0, store->colset_bytes);

  // Copy over bitset, one bit at a time
  for(i = 0; i < fltr->ncols; i++) {
    intocol = file_filter_intocol(fltr, i);
    fromcol = file_filter_fromcol(fltr, i);
    bitset_cpy(store_bitset, intocol, bitset_get(packed_bitset, fromcol));
  }

  // Check colours filtered are actually used in path passed
  for(i = 0, tmp = 0; i < store->colset_bytes; i++) tmp |= store_bitset[i];
  if(tmp == 0) return PATH_NULL;

  // Copy length and path
  const uint8_t *remain_packed = packed + sizeof(PathIndex) + packed_bitset_bytes;
  uint8_t *remain_store = store->next + sizeof(PathIndex) + store->colset_bytes;

  memcpy(remain_store, remain_packed, sizeof(PathLen)+path_nbytes);

  PathIndex match = _path_store_find_or_add_packed(store, last_index, store->next,
                                                   path_nbytes, find, added);

  // Copy path bitset over (may have been found elsewhere)
  uint8_t *match_bitset = store->store + match + sizeof(PathIndex);
  for(i = 0; i < store->colset_bytes; i++) match_bitset[i] |= store_bitset[i];

  return match;
}

// if packed_bases is NULL, uses bases and does packing
// Returns position added to
PathIndex path_store_find_or_add(PathStore *paths, PathIndex last_index,
                                 PathLen len, const Nucleotide *bases,
                                 Orientation orient, Colour colour,
                                 boolean *added)
{
  // We add the path the end of the paths and then check if it is a duplicate
  // this is done because we need to pack the bases into somewhere first anyway
  size_t nbytes = packedpath_len_nbytes(len);
  size_t total_len = packedpath_mem2(paths->colset_bytes, nbytes);
  if(paths->next + total_len >= paths->end) die("Out of memory for paths");

  uint8_t *ptr = paths->next;
  PathLen len_and_orient = packedpath_combine_lenorient(len, orient);

  // write path
  memcpy(ptr, &last_index, sizeof(PathIndex));
  ptr += sizeof(PathIndex);
  memset(ptr, 0, paths->colset_bytes);
  ptr += paths->colset_bytes;
  memcpy(ptr, &len_and_orient, sizeof(PathLen));
  ptr += sizeof(PathLen);
  pack_bases(ptr, bases, len);
  // ptr += nbytes;

  PathIndex match = path_store_find_or_add_packed(paths, last_index, paths->next,
                                                  nbytes, added);

  bitset_set(paths->store+match+sizeof(PathIndex), colour);
  return match;
}

void path_store_fetch_bases(const PathStore *paths, PathIndex pindex,
                            Nucleotide *bases, PathLen len)
{
  uint8_t *ptr = paths->store + pindex + sizeof(PathIndex) +
                 paths->colset_bytes + sizeof(PathLen);
  unpack_bases(ptr, bases, len);
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

  Nucleotide bases[len];
  unpack_bases(seq, bases, len);

  size_t i;
  printf("%8zu: ", (size_t)pindex);
  if(prev == PATH_NULL) printf("    NULL");
  else printf("%8zu", (size_t)prev);
  printf(" (cols:");
  for(i = 0; i < paths->num_of_cols; i++)
    if(bitset_get(colbitset, i)) printf(" %zu", i);
  printf(")[%u]: ", len);
  for(i = 0; i < len; i++)
    putc(dna_nuc_to_char(bases[i]), stdout);
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

  packedpath_get_len_orient(packed-sizeof(PathIndex)-pstore->colset_bytes,
                            pstore->colset_bytes, &plen, &orient);

  Nucleotide bases[plen]; char nucs[plen+1];
  unpack_bases(packed+sizeof(PathLen), bases, plen);
  for(i = 0; i < plen; i++) nucs[i] = dna_nuc_to_char(bases[i]);
  nucs[plen] = '\0';
  status("Path: %zu:%i len %zu %s", (size_t)hkey, orient, (size_t)plen, nucs);
}
