#include "global.h"
#include "db_graph.h"
#include "binary_paths.h"

void path_init(path_t *path)
{
  path->core.prev = PATH_NULL;
  memset(&path->core.colours, 0, sizeof(col_bitset_t));
  path->core.len = 0;
}

void path_alloc(path_t *path)
{
  // printf(" path_allocptr: %zu\n", (size_t)path);
  path_init(path);
  path->bpcap = 16;
  path->cmpcap = 4;
  path->bases = malloc(path->bpcap * sizeof(Nucleotide));
  path->cmpctseq = malloc(path->cmpcap * sizeof(uint8_t));
}

void path_dealloc(path_t *path)
{
  free(path->bases);
  free(path->cmpctseq);
}

void path_to_printf(const path_t *path)
{
  size_t i;
  for(i = 0; i < path->core.len; i++)
    printf(" %c", binary_nuc_to_char(path->bases[i]));
  printf("\n");
}

// {[1:uint64_t prev][N:uint64_t col_bitfield][1:uint32_t len][M:uint8_t data]}..
// prev = PATH_NULL if not set

void binary_paths_init(binary_paths_t *paths, uint8_t *data, size_t size)
{
  binary_paths_t new_paths = {.store = data, .end = data + size,
                              .size = size, .next = data};
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

static inline boolean binary_paths_match(const uint8_t *ptr, const path_t *find)
{
  uint32_t len;
  memcpy(&len, ptr+sizeof(uint64_t)+sizeof(col_bitset_t), sizeof(uint32_t));

  return (len == find->core.len &&
          memcmp(ptr+sizeof(path_core_t), find->cmpctseq,
                 round_bits_to_bytes(len)) == 0);
}

// returns PATH_NULL if not found, otherwise index
static inline uint64_t binary_paths_find(const binary_paths_t *paths,
                                         const path_t *find)
{
  uint64_t last = find->core.prev;
  if(last == PATH_NULL) return false;

  uint8_t *ptr = paths->store + last;
  if(binary_paths_match(ptr, find)) return last;

  while(1)
  {
    memcpy(&last, ptr, sizeof(uint64_t));
    if(last == PATH_NULL) break;

    ptr = paths->store + last;
    if(binary_paths_match(ptr, find)) return last;
  }

  return PATH_NULL;
}

// Returns position added to, or PATH_NULL updated old
uint64_t binary_paths_add(binary_paths_t *paths, const path_t *path, Colour col)
{
  pack_bases(path->cmpctseq, path->bases, path->core.len);
  uint64_t last;

  if(path->core.prev != PATH_NULL &&
     (last = binary_paths_find(paths, path)) != PATH_NULL)
  {
    bitset_set(paths->store+last+sizeof(uint64_t), col);
    return PATH_NULL;
  }

  #ifdef DEBUG
    printf(" Path ADDED\n");
  #endif

  path_to_printf(path);

  size_t len_in_bytes = round_bits_to_bytes(path->core.len*2);
  size_t total_len = sizeof(path_core_t) + len_in_bytes;

  if(paths->next + total_len >= paths->end) die("Out of memory");

  uint8_t *ptr = paths->next;
  uint64_t start = ptr - paths->store;

  // write core
  memcpy(ptr, &path->core, sizeof(path_core_t));
  ptr += sizeof(path_core_t);

  // write bases
  memcpy(ptr, path->cmpctseq, len_in_bytes);
  ptr += len_in_bytes;

  paths->next = ptr;

  return start;
}

void binary_paths_fetch(const binary_paths_t *paths, uint64_t index, path_t *path)
{
  size_t len_in_bytes = round_bits_to_bytes(path->core.len * 2);

  if(len_in_bytes > path->bpcap)
  {
    path->bpcap = ROUNDUP2POW(len_in_bytes);
    path->bases = realloc(path->bases, path->bpcap * sizeof(Nucleotide));
  }

  #ifdef DEBUG
    printf(" Path FETCH\n");
  #endif

  uint8_t *ptr = paths->store + index;

  // read main body
  memcpy(&path->core, ptr, sizeof(path_core_t));
  ptr += sizeof(path_core_t);

  // read bases
  unpack_bases(ptr, path->bases, path->core.len);

  path->pos = 0;
  path->index = index;
}

// Returns 1 on success, 0 otherwise
boolean binary_paths_prev(const binary_paths_t *paths,
                          const path_t *after, path_t *into)
{
  if(after->core.prev == PATH_NULL) return 0;
  binary_paths_fetch(paths, after->core.prev, into);
  return 1;
}
