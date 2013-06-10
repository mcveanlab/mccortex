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
  uint64_t index;
  path_core_t core;
  uint8_t *cmpctseq;
  Nucleotide *bases;
  size_t pos, bpcap, cmpcap;
} path_t;

void path_init(path_t *path);
void path_alloc(path_t *path);
void path_dealloc(path_t *path);

#define path_has_col(path,col) bitset_has((path)->core.colours,col)
#define path_set_col(path,col) bitset_set((path)->core.colours,col)
#define path_del_col(path,col) bitset_del((path)->core.colours,col)

#define PATH_NULL UINT64_MAX

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint32_t len][M:uint8_t data]}..
// N=round_up(NUM_OF_COLOURS/8)
// M=round_up(len/8)

void binary_paths_init(binary_paths_t *paths, uint8_t *data, size_t size);
uint64_t binary_paths_add(binary_paths_t *paths, path_t *path, Colour col);

// unpacks into path_t
void binary_paths_fetch(const binary_paths_t *paths, uint64_t index, path_t *path);

// Returns 1 on success, 0 otherwise
// if found, unpacks into path_t
boolean binary_paths_prev(const binary_paths_t *paths,
                          const path_t *after, path_t *into);

#define binary_paths_del_col(p,idx,col) bitset_del((p)->store+idx+sizeof(uint64_t),col)
#define binary_paths_set_col(p,idx,col) bitset_set((p)->store+idx+sizeof(uint64_t),col)

void binary_paths_dump_path(const path_t *path);

// These are exported for testing only
void pack_bases(uint8_t *ptr, const Nucleotide *bases, size_t len);
void unpack_bases(const uint8_t *ptr, Nucleotide *bases, size_t len);
void check_unpack_capacity(path_t *path, size_t len);
void check_pack_capacity(path_t *path, size_t len);

void binary_paths_dump(const binary_paths_t *paths);

#endif /* BINARY_PATH_H_ */
