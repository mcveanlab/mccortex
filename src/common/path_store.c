#include "global.h"
#include "db_graph.h"
#include "path_store.h"

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}..
// prev = PATH_NULL if not set

void path_store_init(PathStore *paths, uint8_t *data, size_t size,
                       size_t num_of_cols)
{
  uint32_t col_bitset_bytes = round_bits_to_bytes(num_of_cols);

  PathStore new_paths = {.store = data, .end = data + size,
                         .size = size, .next = data,
                         .num_of_cols = num_of_cols,
                         .col_bitset_bytes = col_bitset_bytes,
                         .num_of_paths = 0, .num_kmers_with_paths = 0};

  memcpy(paths, &new_paths, sizeof(PathStore));
}

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


// returns PATH_NULL if not found, otherwise index
// len_in_bytes is length in bytes of bases
static inline PathIndex path_store_find(const PathStore *paths,
                                         PathIndex last_index,
                                         const uint8_t *query,
                                         size_t len_in_bytes)
{
  uint8_t *ptr;
  size_t offset = sizeof(PathIndex) + paths->col_bitset_bytes;
  size_t mem = sizeof(PathLen) + len_in_bytes;

  while(last_index != PATH_NULL)
  {
    ptr = paths->store + last_index;
    if(memcmp(ptr+offset, query+offset, mem) == 0) return last_index;
    memcpy(&last_index, ptr, sizeof(PathIndex));
  }

  return PATH_NULL;
}

// Returns match PathIndex if found, otherwise PATH_NULL
PathIndex validate_appended(PathStore *paths, PathIndex last_index,
                            uint8_t *query, size_t len_in_bytes,
                            boolean *inserted)
{
  PathIndex match;

  if(last_index == PATH_NULL ||
     (match = path_store_find(paths, last_index,
                                query, len_in_bytes)) == PATH_NULL)
  {
    // Accept new path
    size_t mem = sizeof(PathIndex) + paths->col_bitset_bytes +
                 sizeof(PathLen) + len_in_bytes;

    if(query != paths->next)
      memcpy(paths->next, query, mem);

    match = paths->next - paths->store;
    paths->next += mem;
    paths->num_of_paths++;
    paths->num_kmers_with_paths += (last_index == PATH_NULL);
    *inserted = true;
  }
  else *inserted = false;

  return match;
}

PathIndex path_store_add2(PathStore *paths, PathIndex last_index,
                            uint8_t *packed)
{
  PathLen len;
  boolean inserted = false;
  PathIndex idx;
  size_t i;

  memcpy(&len, packed+sizeof(PathIndex)+paths->col_bitset_bytes, sizeof(PathLen));
  size_t nbytes = (len+3)/4;
  idx = validate_appended(paths, last_index, packed, nbytes, &inserted);
  if(!inserted) {
    uint8_t *dst = paths->store + idx + sizeof(PathIndex);
    uint8_t *src = packed + sizeof(PathIndex);
    for(i = 0; i < paths->col_bitset_bytes; i++) dst[i] |= src[i];
  }
  return inserted ? idx : PATH_NULL;
}

// if packed_bases is NULL, uses bases and does packing
// Returns position added to
PathIndex path_store_add(PathStore *paths, PathIndex last_index,
                           PathLen len, const Nucleotide *bases,
                           Orientation orient, Colour colour)
{
  // We add the path the end of the paths and then check if it is a duplicate
  // this is done because we need to pack the bases into somewhere first anyway
  size_t nbytes = (len+3)/4;
  size_t total_len = path_mem(paths->col_bitset_bytes, len);

  if(paths->next + total_len >= paths->end) die("Out of memory");

  uint8_t *ptr = paths->next;

  uint8_t colbitset[paths->col_bitset_bytes];
  memset(colbitset, 0, paths->col_bitset_bytes);

  PathLen len_and_orient = len | (orient << PATH_LEN_BITS);

  // write path
  memcpy(ptr, &last_index, sizeof(PathIndex));
  ptr += sizeof(PathIndex);
  memcpy(ptr, colbitset, paths->col_bitset_bytes);
  ptr += paths->col_bitset_bytes;
  memcpy(ptr, &len_and_orient, sizeof(PathLen));
  ptr += sizeof(PathLen);
  pack_bases(ptr, bases, len);
  // ptr += nbytes;

  boolean inserted;
  PathIndex match;
  match = validate_appended(paths, last_index, paths->next, nbytes, &inserted);
  bitset_set(paths->store+match+sizeof(PathIndex), colour);

  return inserted ? match : PATH_NULL;
}

PathIndex path_store_prev(const PathStore *paths, PathIndex index)
{
  PathIndex prev;
  memcpy(&prev, paths->store + index, sizeof(PathIndex));
  return prev;
}

void path_store_len_orient(const PathStore *paths, PathIndex index,
                             PathLen *len, Orientation *orient)
{
  PathLen d;
  memcpy(&d, paths->store + index + sizeof(PathIndex) + paths->col_bitset_bytes,
         sizeof(PathLen));
  *len = d & ~(1<<PATH_LEN_BITS);
  *orient = d >> PATH_LEN_BITS;
}

void path_store_fetch(const PathStore *paths, PathIndex index,
                        Nucleotide *bases, PathLen len)
{
  uint8_t *ptr = paths->store + index + sizeof(PathIndex) +
                 paths->col_bitset_bytes + sizeof(PathLen);
  unpack_bases(ptr, bases, len);
}

size_t path_store_size(const PathStore *paths, PathIndex index)
{
  PathLen len;
  Orientation orient;
  path_store_len_orient(paths, index, &len, &orient);
  return sizeof(PathIndex) + paths->col_bitset_bytes + sizeof(PathLen) +
         (len+3)/4;
}

void path_store_print_path(const PathStore *paths, PathIndex index)
{
  PathIndex prev;
  PathLen len;
  Orientation orient;

  prev = path_store_prev(paths, index);
  path_store_len_orient(paths, index, &len, &orient);
  uint8_t *colbitset = paths->store + index + sizeof(PathIndex);
  uint8_t *data = colbitset + paths->col_bitset_bytes + sizeof(PathLen);

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
  PathIndex index, store_size = paths->next - paths->store;

  for(index = 0; index < store_size; index += path_store_size(paths, index))
    path_store_print_path(paths, index);
}
