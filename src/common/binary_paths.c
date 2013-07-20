#include "global.h"
#include "db_graph.h"
#include "binary_paths.h"

void path_init(path_t *path, uint32_t num_of_cols)
{
  path_t tmp = {.num_of_cols = num_of_cols,
                .col_bitset_bytes = round_bits_to_bytes(num_of_cols),
                .index = PATH_NULL, .prev = PATH_NULL, .len = 0,
                .data = path->data, .datacap = path->datacap,
                .bases = path->bases, .bpcap = path->bpcap};
  memcpy(path, &tmp, sizeof(path_t));
  // Clear colour bitset
  memset(path->data, 0, path->col_bitset_bytes);
}

void path_alloc(path_t *path, uint32_t num_of_cols)
{
  path->bpcap = path->datacap = 16;
  path->bases = malloc(path->bpcap * sizeof(Nucleotide));
  path->data = malloc(path->datacap * sizeof(uint8_t));
  path_init(path, num_of_cols);
}

void path_dealloc(path_t *path)
{
  free(path->data);
  free(path->bases);
}

void path_to_printf(const path_t *path)
{
  size_t i;
  for(i = 0; i < path->len; i++)
    printf(" %c", binary_nuc_to_char(path->bases[i]));
  printf("\n");
}

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint32_t len][M:uint8_t data]}..
// prev = PATH_NULL if not set

void binary_paths_init(binary_paths_t *paths, uint8_t *data, size_t size,
                       size_t num_of_cols)
{
  binary_paths_t new_paths = {.store = data, .end = data + size,
                              .size = size, .next = data,
                              .num_of_paths = 0, .num_of_cols = num_of_cols};
  memcpy(paths, &new_paths, sizeof(binary_paths_t));
}

void pack_bases(uint8_t *ptr, const Nucleotide *bases, size_t len)
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

void unpack_bases(const uint8_t *ptr, Nucleotide *bases, size_t len)
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

void check_unpack_capacity(path_t *path, size_t len)
{
  if(len > path->bpcap)
  {
    path->bpcap = ROUNDUP2POW(len);
    path->bases = realloc(path->bases, path->bpcap * sizeof(Nucleotide));
  }
}

void check_pack_capacity(path_t *path, size_t len)
{
  size_t len_in_bytes = round_bits_to_bytes(len*2) + path->col_bitset_bytes;

  if(len_in_bytes > path->datacap)
  {
    path->datacap = ROUNDUP2POW(len_in_bytes);
    path->data = realloc(path->data, path->datacap * sizeof(uint8_t));
  }
}

static inline boolean binary_paths_match(const uint8_t *ptr, const path_t *find)
{
  size_t offset = sizeof(uint64_t) + find->col_bitset_bytes;
  size_t len_in_bytes = round_bits_to_bytes(find->len*2);
  return (memcmp(ptr+offset, find->data, len_in_bytes) == 0);
}

// returns PATH_NULL if not found, otherwise index
static inline uint64_t binary_paths_find(const binary_paths_t *paths,
                                         const path_t *find)
{
  uint8_t *ptr;
  uint64_t last = find->prev;
  if(last == PATH_NULL) return PATH_NULL;

  while(last != PATH_NULL)
  {
    ptr = paths->store + last;
    if(binary_paths_match(ptr, find)) return last;
    memcpy(&last, ptr, sizeof(uint64_t));
  }

  return PATH_NULL;
}

// Returns position added to, or PATH_NULL if updated an old entry
uint64_t binary_paths_add(binary_paths_t *paths, path_t *path, Colour colour)
{
  size_t len_in_bytes = round_bits_to_bytes(path->len*2);

  // Check capacity
  check_pack_capacity(path, path->len);

  pack_bases(path_packed_bases(path), path->bases, path->len);
  uint64_t last;

  if(path->prev != PATH_NULL &&
     (last = binary_paths_find(paths, path)) != PATH_NULL)
  {
    bitset_set(paths->store+last+sizeof(uint64_t), colour);
    return PATH_NULL;
  }

  size_t total_len = sizeof(uint64_t) + path->col_bitset_bytes +
                     sizeof(uint32_t) + len_in_bytes;

  if(paths->next + total_len >= paths->end) die("Out of memory");

  uint8_t *ptr = paths->next;
  uint64_t start = ptr - paths->store;
  path->index = start;

  #ifdef DEBUG
    binary_paths_dump_path(path);
  #endif

  // write path
  memcpy(ptr, &path->prev, sizeof(uint64_t));
  ptr += sizeof(uint64_t);
  memcpy(ptr, &path->data, path->col_bitset_bytes);
  ptr += path->col_bitset_bytes;
  memcpy(ptr, &path->len, sizeof(uint32_t));
  ptr += sizeof(uint32_t);
  memcpy(ptr, path_packed_bases(path), len_in_bytes);
  ptr += len_in_bytes;

  paths->next = ptr;
  paths->num_of_paths++;

  return start;
}

// this function doesn't copy compacted bases to path->data but it does unpack
// into path->bases
void binary_paths_fetch(const binary_paths_t *paths, uint64_t index, path_t *path)
{
  uint8_t *ptr = paths->store + index;
  path->index = index;
  memcpy(&path->prev, ptr, sizeof(uint64_t));
  ptr += sizeof(uint64_t);
  memcpy(path->data, ptr, path->col_bitset_bytes);
  ptr += path->col_bitset_bytes;
  memcpy(&path->len, ptr, sizeof(uint32_t));
  ptr += sizeof(uint32_t);
  check_unpack_capacity(path, path->len);
  unpack_bases(ptr, path->bases, path->len);
  path->pos = 0;
}

// Returns 1 on success, 0 otherwise
boolean binary_paths_prev(const binary_paths_t *paths,
                          const path_t *after, path_t *into)
{
  if(after->prev == PATH_NULL) return 0;
  binary_paths_fetch(paths, after->prev, into);
  return 1;
}

void binary_paths_dump_path(const path_t *path)
{
  size_t i;
  printf("%8zu: ", (size_t)path->index);
  if(path->prev == PATH_NULL) printf("    NULL");
  else printf("%8zu", (size_t)path->prev);
  printf(" (cols:");
  for(i = 0; i < path->num_of_cols; i++)
    if(bitset_has(path->data, i)) printf(" %zu", i);
  printf(")[%zu/%u]: ", path->pos, path->len);
  for(i = 0; i < path->len; i++)
    putc(binary_nuc_to_char(path->bases[i]), stdout);
  putc('\n', stdout);
}

void binary_paths_dump(const binary_paths_t *paths)
{
  path_t tmp;
  uint64_t index;
  path_alloc(&tmp, paths->num_of_cols);
  const uint64_t len = paths->next - paths->store;

  for(index = 0; index < len; index += path_size(&tmp))
  {
    binary_paths_fetch(paths, index, &tmp);
    binary_paths_dump_path(&tmp);
  }

  path_dealloc(&tmp);
}
