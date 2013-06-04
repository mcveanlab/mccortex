#ifndef BINARY_PATH_H_
#define BINARY_PATH_H_

#include "binary_kmer.h"

typedef uint8_t col_bitset_t[round_bits_to_bytes(NUM_OF_COLOURS)];

typedef struct path_core_t path_core_t;

struct path_core_t
{
  uint64_t prev;
  col_bitset_t colours;
  uint32_t len;
} __attribute__((packed));

typedef struct
{
  path_core_t core;
  uint8_t *cmpctseq;
  Nucleotide *bases;
  size_t pos, bpcap, cmpcap;
} path_t;

void path_init(path_t *path);
void path_alloc(path_t *path);
void path_dealloc(path_t *path);

#define path_has_col(path,col) bitset_has((path)->colours,col)
#define path_set_col(path,col) bitset_set((path)->colours,col)
#define path_del_col(path,col) bitset_del((path)->colours,col)

typedef struct {
  uint8_t *const store, *const end;
  const size_t size;
  uint8_t *next;
} binary_paths_t;

#define PATH_NULL UINT64_MAX

// {[1:uint64_t prev][N:uint64_t col_bitfield][1:uint32_t len][M:uint8_t data]}..

void binary_paths_init(binary_paths_t *paths, uint8_t *data, size_t size);
uint64_t binary_paths_add(binary_paths_t *paths, const path_t *path, Colour col);
void binary_paths_fetch(const binary_paths_t *paths, uint64_t index, path_t *path);

// Returns 1 on success, 0 otherwise
boolean binary_paths_prev(const binary_paths_t *paths,
                          const path_t *after, path_t *into);

void pack_bases(uint8_t *ptr, const Nucleotide *bases, size_t len);
void unpack_bases(const uint8_t *ptr, Nucleotide *bases, size_t len);

#endif /* BINARY_PATH_H_ */
