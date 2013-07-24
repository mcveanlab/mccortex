#ifndef BINARY_PATH_H_
#define BINARY_PATH_H_

#include "binary_kmer.h"

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint32_t len][M:uint8_t data]}..
// N=round_up(num_of_colours/8)
// M=round_up(len/4)

void binary_paths_init(PathStore *paths, uint8_t *data, size_t size,
                       size_t num_of_cols);

// Add
PathIndex binary_paths_add2(PathStore *paths, PathIndex last_index,
                            uint8_t *packed);

PathIndex binary_paths_add(PathStore *paths, PathIndex last_index,
                           PathLen len, const Nucleotide *bases,
                           Orientation orient, Colour colour);

// Fetch
PathIndex binary_paths_prev(const PathStore *paths, PathIndex index);
void binary_paths_len_orient(const PathStore *paths, PathIndex index,
                             PathLen *len, Orientation *orient);
void binary_paths_fetch(const PathStore *paths, PathIndex index,
                        Nucleotide *bases, PathLen len);


#define binary_paths_del_col(p,idx,col) bitset_del((p)->store+idx+sizeof(PathIndex),col)
#define binary_paths_set_col(p,idx,col) bitset_set((p)->store+idx+sizeof(PathIndex),col)
#define binary_paths_has_col(p,idx,col) bitset_has((p)->store+idx+sizeof(PathIndex),col)

// Print
void binary_paths_dump_path(const PathStore *paths, PathIndex index);
void binary_paths_dump(const PathStore *paths);

#endif /* BINARY_PATH_H_ */
