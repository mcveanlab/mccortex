#ifndef BINARY_PATH_H_
#define BINARY_PATH_H_

#include "binary_kmer.h"

typedef struct
{
  const uint32_t num_of_cols, col_bitset_bytes;
  // data is <col_bitset><bases>
  uint8_t *data;
  uint64_t index, prev;
  uint32_t len;
  Nucleotide *bases;
  size_t pos, bpcap, datacap;
} path_t;

#define PATH_NULL UINT64_MAX

#define path_packed_bases(path) ((path)->data+(path)->col_bitset_bytes)

void path_init(path_t *path, uint32_t num_of_cols);
void path_alloc(path_t *path, uint32_t num_of_cols);
void path_dealloc(path_t *path);

#define path_has_col(path,col) bitset_has((path)->data,col)
#define path_set_col(path,col) bitset_set((path)->data,col)
#define path_del_col(path,col) bitset_del((path)->data,col)

#define path_size(path) (8+(path)->col_bitset_bytes+4+round_bits_to_bytes((path)->len*2))

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint32_t len][M:uint8_t data]}..
// N=round_up(num_of_colours/8)
// M=round_up(len/8)

void binary_paths_init(binary_paths_t *paths, uint8_t *data, size_t size,
                       size_t num_of_cols);
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
