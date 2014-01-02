#ifndef BINARY_PATH_H_
#define BINARY_PATH_H_

#include "binary_kmer.h"
#include "file_filter.h"

// Types PathIndex, PathLen are defined in graph_typedef.h
// typedef uint64_t PathIndex;
// typedef uint16_t PathLen;
#define PATH_NULL UINT64_MAX
#define PATH_LEN_BITS 15
#define MAX_PATHLEN ((1U<<PATH_LEN_BITS)-1)

// Data structure:
//   {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}..
// N=round_up(num_of_colours/8), M=round_up(len/4)

// Initialise the PathStore
void path_store_init(PathStore *paths, uint8_t *data, size_t size, size_t ncols);
void path_store_resize(PathStore *paths, size_t size);

// Find or add a path into the PathStore
// last_index is index of the last path belonging to the kmer which owns the
// new path that is to be inserted
// Returns match PathIndex if found, otherwise PATH_NULL
PathIndex path_store_find_or_add_packed(PathStore *paths, PathIndex last_index,
                                        const uint8_t *packed, size_t path_nbytes,
                                        boolean *inserted);

// Add a PackedPath, using a FileFilter to reduce to a subset of colours
// `find` Specifies if we should try to find a duplicate first
// Returns PATH_NULL if no colours set in colour subset
PathIndex path_store_find_or_add_packed2(PathStore *store, PathIndex last_index,
                                         const uint8_t *packed, size_t path_nbytes,
                                         const FileFilter *fltr, boolean find,
                                         boolean *added);

// Add to PathStore
PathIndex path_store_find_or_add(PathStore *paths, PathIndex last_index,
                                 PathLen len, const Nucleotide *bases,
                                 Orientation orient, Colour colour,
                                 boolean *added);

// Print
void path_store_print_path(const PathStore *paths, PathIndex index);
void path_store_print_all(const PathStore *paths);

// Fetch sequence
void path_store_fetch_bases(const PathStore *paths, PathIndex index,
                            Nucleotide *bases, PathLen len);

// If compatible, a FileFilter can be read straight into a PathStore without
// parsing each path, one-by-one (much faster!)
#define path_store_fltr_compatible(st,fltr) \
        ((fltr)->nofilter && \
         roundup_bits2bytes((fltr)->filencols) == (st)->colset_bytes)

//
// Functions for PackedPaths
//
// PackedPath doesn't have a struct, instead is laid out as:
//   {[1:uint64_t prev][N:uint8_t col_bitfield][1:uint16_t len][M:uint8_t data]}
// Where N=round_up(num_of_colours/8), M=round_up(len/4)
// len holds both the length (lower 15 bits), and orientations (top bit)

// Number of bytes needed to store n bases (2 bits per base, 4 per byte)
#define packedpath_len_nbytes(nbases) (((nbases)+3)/4)

// Number of bytes a packed path takes up
#define packedpath_mem2(colbytes,nbytes) \
        (sizeof(PathIndex)+(colbytes)+sizeof(PathLen)+(nbytes))
#define packedpath_mem(p,colbytes) \
        packedpath_mem2(colbytes,packedpath_pbytes(p,colbytes))

#define packedpath_lenword2(w) ((w) & ~(1UL<<PATH_LEN_BITS))
#define packedpath_orword2(w) ((w) >> PATH_LEN_BITS)

#define packedpath_lenword(p,colbytes) (*(PathLen*)((p)+sizeof(PathIndex)+colbytes))
#define packedpath_len(p,colbytes) packedpath_lenword2(packedpath_lenword(p,colbytes))
#define packedpath_orient(p,colbytes) packedpath_orword2(packedpath_lenword(p,colbytes))
#define packedpath_pbytes(p,colbytes) packedpath_len_nbytes(packedpath_len(p,colbytes))
#define packedpath_prev(p) (*(PathIndex*)(p))
#define packedpath_colset(p) ((p)+sizeof(PathIndex))
#define packedpath_path(p,colbits) ((p)+sizeof(PathIndex)+colbits+sizeof(PathLen))

#define packedpath_del_col(p,idx,col) bitset_del((p)+sizeof(PathIndex),col)
#define packedpath_set_col(p,idx,col) bitset_set((p)+sizeof(PathIndex),col)
#define packedpath_has_col(p,idx,col) bitset_get((p)+sizeof(PathIndex),col)

static inline PathLen packedpack_len_orient(const uint8_t *packed,
                                            const PathStore *store,
                                            PathLen *len, Orientation *orient)
{
  PathLen d = packedpath_lenword(packed, store->colset_bytes);
  *len = packedpath_lenword2(d);
  *orient = packedpath_orword2(d);
  return d;
}

#endif /* BINARY_PATH_H_ */
