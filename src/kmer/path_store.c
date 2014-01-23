#include "global.h"
#include "util.h"
#include "db_graph.h"
#include "path_store.h"

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}..
// prev = PATH_NULL if not set

void path_store_init(PathStore *paths, uint8_t *data, size_t size,
                       size_t num_of_cols)
{
  size_t colset_bytes = roundup_bits2bytes(num_of_cols);

  char memstr[100]; bytes_to_str(size, 1, memstr);
  status("[paths] Setting up path store to use %s", memstr);

  PathStore new_paths = {.store = data, .end = data + size,
                         .size = size, .next = data,
                         .num_of_cols = num_of_cols,
                         .colset_bytes = colset_bytes,
                         .num_of_paths = 0, .num_kmers_with_paths = 0};

  memcpy(paths, &new_paths, sizeof(PathStore));
}

void path_store_resize(PathStore *paths, size_t size)
{
  char memstr[100]; bytes_to_str(size, 1, memstr);
  status("[paths] Resizing path store to use %s", memstr);

  PathStore new_paths = {.store = paths->store, .end = paths->store + size,
                         .size = size, .next = paths->next,
                         .num_of_cols = paths->num_of_cols,
                         .colset_bytes = paths->colset_bytes,
                         .num_of_paths = paths->num_of_paths,
                         .num_kmers_with_paths = paths->num_kmers_with_paths};

  memcpy(paths, &new_paths, sizeof(PathStore));
}

// Convert from unpacked representation (1 base per byte) to packed
// representation (4 bases per byte)
void pack_bases(uint8_t *ptr, const Nucleotide *bases, size_t len)
{
  size_t full_bytes = len/4;
  const uint8_t *endptr = ptr+full_bytes;

  for(; ptr < endptr; ptr++)
    *ptr = bases[0] | (bases[1]<<2) | (bases[2]<<4) | (bases[3]<<6);

  // Do last byte
  if(len & 3) *ptr = 0;
  switch(len & 3) {
    case 3: *ptr = bases[--len];
    case 2: *ptr = (*ptr<<2) | bases[--len];
    case 1: *ptr = (*ptr<<2) | bases[--len];
  }
}

// Convert from compact representation (4 bases per byte) to unpacked
// representation (1 base per byte)
void unpack_bases(const uint8_t *ptr, Nucleotide *bases, size_t len)
{
  size_t i, full_bytes = len/4;
  const uint8_t *endptr = ptr+full_bytes;

  for(i = 0; ptr < endptr; ptr++) {
    bases[i++] =  (*ptr)     & 3;
    bases[i++] = ((*ptr)>>2) & 3;
    bases[i++] = ((*ptr)>>4) & 3;
    bases[i++] = ((*ptr)>>6);
  }

  // Do last byte
  switch(len & 3) {
    case 3: bases[i++] = (*ptr>>4) & 3;
    case 2: bases[i++] = (*ptr>>2) & 3;
    case 1: bases[i++] = (*ptr)    & 3;
  }
}

// Copy a packed path from one place in memory to another, applying left shift
// Shifting by N bases results in N fewer bases in output
void packed_cpy(uint8_t *restrict dst, const uint8_t *restrict src,
                size_t shift, size_t len_bases)
{
  assert(shift < 4);
  size_t i, nbytes = (len_bases+3)/4;

  if(shift >= len_bases) return;
  if(!shift) { memcpy(dst, src, nbytes); return; }

  switch(shift) {
    case 3: for(i=0;i+1<nbytes;i++) dst[i] = (src[i]>>6) | (src[i+1]<<2); break;
    case 2: for(i=0;i+1<nbytes;i++) dst[i] = (src[i]>>4) | (src[i+1]<<4); break;
    case 1: for(i=0;i+1<nbytes;i++) dst[i] = (src[i]>>2) | (src[i+1]<<6); break;
  }
  dst[nbytes-1] = src[nbytes-1] >> shift;
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

  PathIndex index = (PathIndex)(store->next - store->store);
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
  PathLen len_and_orient = len | (PathLen)(orient << PATH_LEN_BITS);

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
    if(bitset_get(colbitset, i)) printf(" %zu", i);
  printf(")[%u]: ", len);
  for(i = 0; i < len; i++)
    putc(dna_nuc_to_char(bases[i]), stdout);
  putc('\n', stdout);
}

void path_store_print_all(const PathStore *paths)
{
  PathIndex index = 0, store_size = (PathIndex)(paths->next - paths->store);

  while(index < store_size) {
    path_store_print_path(paths, index);
    index += packedpath_mem(paths->store+index, paths->colset_bytes);
  }
}
