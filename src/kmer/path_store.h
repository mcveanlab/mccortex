#ifndef BINARY_PATH_H_
#define BINARY_PATH_H_

#include "binary_kmer.h"
#include "file_filter.h"
#include "cortex_types.h"
#include "packed_path.h"

// Extra padding to avoid reading bad memory
#define PSTORE_PADDING 16

typedef struct {
  uint8_t *const store, *const end;
  const size_t size;
  const size_t num_of_cols, colset_bytes;
  uint8_t *next;
  size_t num_of_paths, num_kmers_with_paths;
  // Temporary data used for merging
  uint8_t *const tmpdata;
  const size_t tmpsize;
} PathStore;

// Initialise the PathStore
void path_store_alloc(PathStore *paths, size_t size, size_t tmpsize, size_t ncols);
// Once tmp has been used for merging, it can be reclaimed to use generally
void path_store_reclaim_tmp(PathStore *paths);

// Release memory
void path_store_dealloc(PathStore *paths);

// Find a path
// returns PATH_NULL if not found, otherwise index
// path_nbytes is length in bytes of bases = (num bases + 3)/4
// query is <PathLen><PackedSeq>
PathIndex path_store_find(const PathStore *paths, PathIndex last_index,
                          const uint8_t *query, size_t path_nbytes);

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

// Fetch sequence
void path_store_fetch_bases(const PathStore *paths, PathIndex index,
                            Nucleotide *bases, PathLen len);

// If compatible, a FileFilter can be read straight into a PathStore without
// parsing each path, one-by-one (much faster!)
#define path_store_fltr_compatible(st,fltr) \
        ((fltr)->nofilter && \
         roundup_bits2bytes((fltr)->filencols) == (st)->colset_bytes)

// Print
void path_store_print_path(const PathStore *paths, PathIndex index);
void path_store_print_all(const PathStore *paths);

// packed points to <PathLen><PackedSeq>
void print_path(hkey_t hkey, const uint8_t *packed, const PathStore *pstore);

#endif /* BINARY_PATH_H_ */
