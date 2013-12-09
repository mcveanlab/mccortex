#include "global.h"
#include "util.h"
#include "db_graph.h"
#include "path_store.h"

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}..
// prev = PATH_NULL if not set

void path_store_init(PathStore *paths, uint8_t *data, size_t size,
                       size_t num_of_cols)
{
  uint32_t colset_bytes = round_bits_to_bytes(num_of_cols);

  PathStore new_paths = {.store = data, .end = data + size,
                         .size = size, .next = data,
                         .num_of_cols = num_of_cols,
                         .colset_bytes = colset_bytes,
                         .num_of_paths = 0, .num_kmers_with_paths = 0};

  memcpy(paths, &new_paths, sizeof(PathStore));
}

void path_store_resize(PathStore *paths, size_t size)
{
  PathStore new_paths = {.store = paths->store, .end = paths->store + size,
                         .size = size, .next = paths->next,
                         .num_of_cols = paths->num_of_cols,
                         .colset_bytes = paths->colset_bytes,
                         .num_of_paths = paths->num_of_paths,
                         .num_kmers_with_paths = paths->num_kmers_with_paths};

  memcpy(paths, &new_paths, sizeof(PathStore));
}

// Convert from unpacked representation (1 bas per byte) to packed
// representation (4 bases per byte)
static inline void pack_bases(uint8_t *ptr, const Nucleotide *bases, size_t len)
{
  uint8_t tmp;
  size_t i, j, k = 0, full_bytes = len/4;
  for(i = 0; i < full_bytes; i++)
  {
    tmp = 0; k += 4;
    for(j = 0; j < 4; j++) {
      tmp <<= 2;
      tmp |= bases[--k];
    }
    k += 4;
    *ptr = tmp;
    ptr++;
  }

  // Do last byte
  size_t quotient = full_bytes*4;
  if(quotient < len)
  {
    tmp = 0; k = len;
    while(k > quotient) {
      tmp <<= 2;
      tmp |= bases[--k];
    }
    *ptr = tmp;
  }
}

// Convert from compact representation (4 bases per byte) to unpacked
// representation (1 base per byte)
static inline void unpack_bases(const uint8_t *ptr, Nucleotide *bases, size_t len)
{
  uint8_t tmp;
  size_t i, j, k = 0, full_bytes = len/4;
  for(i = 0; i < full_bytes; i++)
  {
    tmp = *ptr;
    for(j = 0; j < 4; j++) {
      bases[k++] = tmp & 0x3;
      tmp >>= 2;
    }
    ptr++;
  }

  // Do last byte
  tmp = *ptr;
  while(k < len) {
    bases[k++] = tmp & 0x3;
    tmp >>= 2;
  }
}

// Find a path
// returns PATH_NULL if not found, otherwise index
// path_nbytes is length in bytes of bases = (num bases + 3)/4
static inline PathIndex path_store_find(const PathStore *paths,
                                        PathIndex last_index,
                                        const uint8_t *query,
                                        size_t path_nbytes)
{
  uint8_t *packed;
  size_t offset = sizeof(PathIndex) + paths->colset_bytes;
  size_t mem = sizeof(PathLen) + path_nbytes;

  while(last_index != PATH_NULL)
  {
    packed = paths->store + last_index;
    if(memcmp(packed+offset, query+offset, mem) == 0) return last_index;
    last_index = packedpath_prev(packed);
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

  PathIndex index = store->next - store->store;
  store->next += mem;
  store->num_of_paths++;
  store->num_kmers_with_paths += (last_index == PATH_NULL);
  return index;
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
  PathIndex match = path_store_find(paths, last_index, packed, path_nbytes);

  if(match == PATH_NULL)
  {
    match = path_store_add_packed(paths, last_index, packed, path_nbytes);
    *inserted = true;
  }
  else {
    // Already in path store, just update colour bitset
    uint8_t *dst = packedpath_colset(paths->store + match);
    const uint8_t *src = packedpath_colset(packed);
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
  size_t packed_bitset_bytes = round_bits_to_bytes(fltr->filencols);

  if(path_store_fltr_compatible(store,fltr)) {
    return _path_store_find_or_add_packed(store, last_index, packed,
                                         path_nbytes, find, added);
  }

  // Check we have enough memory to add
  size_t mem = packedpath_mem2(store->colset_bytes, path_nbytes);
  if(store->next + mem >= store->end) die("Out of memory for paths");

  uint8_t *ptr = store->next;
  memcpy(ptr, &last_index, sizeof(PathIndex));
  ptr += sizeof(PathIndex);

  // Clear memory for colour bitset
  memset(ptr, 0, store->colset_bytes);

  // Copy over bitset, one bit at a time
  const uint8_t *packed_bitset = packed + sizeof(PathIndex);
  size_t i, tocol, fromcol;
  for(i = 0; i < fltr->ncols; i++) {
    tocol = file_filter_intocol(fltr, i);
    fromcol = fltr->cols[i];
    bitset_cpy(ptr, tocol, bitset_has(packed_bitset, fromcol));
  }

  // Check colours are used
  uint8_t tmp = 0;
  for(i = 0; i < store->colset_bytes; i++) tmp |= ptr[i];
  if(tmp == 0) return PATH_NULL;

  // Copy length and path
  ptr += store->colset_bytes;
  packed += sizeof(PathIndex) + packed_bitset_bytes;
  memcpy(ptr, packed, sizeof(PathLen)+path_nbytes);

  return _path_store_find_or_add_packed(store, last_index, store->next,
                                        path_nbytes, find, added);
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
  PathLen len_and_orient = len | (orient << PATH_LEN_BITS);

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

void path_store_fetch_bases(const PathStore *paths, PathIndex index,
                            Nucleotide *bases, PathLen len)
{
  uint8_t *ptr = paths->store + index + sizeof(PathIndex) +
                 paths->colset_bytes + sizeof(PathLen);
  unpack_bases(ptr, bases, len);
}

void path_store_print_path(const PathStore *paths, PathIndex index)
{
  PathIndex prev;
  PathLen len;
  Orientation orient;

  prev = packedpath_prev(paths->store + index);
  packedpack_len_orient(paths->store+index, paths, &len, &orient);

  const uint8_t *packed = paths->store + index;
  const uint8_t *colbitset = packedpath_colset(packed);
  const uint8_t *data = packedpath_path(packed, paths->colset_bytes);

  Nucleotide bases[len];
  unpack_bases(data, bases, len);

  size_t i;
  printf("%8zu: ", (size_t)index);
  if(prev == PATH_NULL) printf("    NULL");
  else printf("%8zu", (size_t)prev);
  printf(" (cols:");
  for(i = 0; i < paths->num_of_cols; i++)
    if(bitset_has(colbitset, i)) printf(" %zu", i);
  printf(")[%u]: ", len);
  for(i = 0; i < len; i++)
    putc(binary_nuc_to_char(bases[i]), stdout);
  putc('\n', stdout);
}

void path_store_print_all(const PathStore *paths)
{
  PathIndex index = 0, store_size = paths->next - paths->store;

  while(index < store_size) {
    path_store_print_path(paths, index);
    index += packedpath_mem(paths->store+index, paths->colset_bytes);
  }
}
