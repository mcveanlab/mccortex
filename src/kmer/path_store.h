#ifndef BINARY_PATH_H_
#define BINARY_PATH_H_

#include "binary_kmer.h"

// {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint32_t len][M:uint8_t data]}..
// N=round_up(num_of_colours/8)
// M=round_up(len/4)

void path_store_init(PathStore *paths, uint8_t *data, size_t size, size_t ncols);

// Add
PathIndex path_store_add2(PathStore *paths, PathIndex last_index,
                            uint8_t *packed);

PathIndex path_store_add(PathStore *paths, PathIndex last_index,
                           PathLen len, const Nucleotide *bases,
                           Orientation orient, Colour colour);

// Fetch
PathIndex path_store_prev(const PathStore *paths, PathIndex index);
void path_store_len_orient(const PathStore *paths, PathIndex index,
                             PathLen *len, Orientation *orient);
void path_store_fetch(const PathStore *paths, PathIndex index,
                        Nucleotide *bases, PathLen len);

#define path_len_in_bytes(nbases) (((nbases)+3)/4)
#define path_mem(colbytes,nbases) \
        (sizeof(PathIndex)+colbytes+sizeof(PathLen)+((nbases)+3)/4)

#define path_store_del_col(p,idx,col) bitset_del((p)->store+idx+sizeof(PathIndex),col)
#define path_store_set_col(p,idx,col) bitset_set((p)->store+idx+sizeof(PathIndex),col)
#define path_store_has_col(p,idx,col) bitset_has((p)->store+idx+sizeof(PathIndex),col)

// Print
void path_store_print_path(const PathStore *paths, PathIndex index);
void path_store_print_all(const PathStore *paths);

#endif /* BINARY_PATH_H_ */
