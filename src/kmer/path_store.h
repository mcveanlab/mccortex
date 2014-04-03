#ifndef BINARY_PATH_H_
#define BINARY_PATH_H_

#include "binary_kmer.h"
#include "file_filter.h"
#include "cortex_types.h"
#include "packed_path.h"
#include "path_hash.h"

// Extra padding to avoid reading bad memory
#define PSTORE_PADDING 16

typedef struct
{
  // Constants for this instance
  uint8_t *const store, *end, *next;
  const size_t num_of_cols, colset_bytes;

  size_t num_of_paths, num_kmers_with_paths, num_col_paths;

  // Temporary data used for merging
  uint8_t *tmpstore;
  size_t tmpsize;

  // Kmer pointers
  // kmer_paths is NULL if we don't want GraphWalker to use paths
  // kmer_paths_update may point to kmer_paths, elsewhere or NULL
  //
  // kmer_paths is used in loading paths and traversing the graph
  // kmer_paths_update is used for generate_paths.c and dumping ctp files
  PathIndex *kmer_paths, *kmer_paths_update;

  // Multithreaded writing
  uint8_t *kmer_locks;

  PathHash phash;
} PathStore;

//
// Pointer to linked list for each kmer
//
#define pstore_paths(pstore,node) \
        ((pstore)->kmer_paths[(node)])
#define pstore_paths_updated(pstore,node) \
        ((pstore)->kmer_paths_update[(node)])
#define pstore_paths_update_volptr(pstore,node) \
        ((volatile PathIndex *)&(pstore)->kmer_paths_update[(node)])

// Initialise the PathStore
// use_path_hash must be true if you are adding paths manual
//   (i.e. adding paths anyway other than loading from a file)
void path_store_alloc(PathStore *ps, size_t mem, bool use_path_hash,
                      size_t kmers_in_hash, size_t ncols);

// Release memory
void path_store_dealloc(PathStore *paths);


// Set up temporary memory for merging PathStores
void path_store_setup_tmp(PathStore *ps, size_t tmp_mem);

// Once tmp has been used for merging, it can be reclaimed to use generally
void path_store_release_tmp(PathStore *ps);

// Free ps->kmer_paths, set kmer_paths to point to kmer_paths_update.
void path_store_combine_updated_paths(PathStore *pstore);

// Find a path
// returns PATH_NULL if not found, otherwise index
// path_nbytes is length in bytes of bases = (num bases + 3)/4
// query is <PathLen><PackedSeq>
PathIndex path_store_find(const PathStore *paths, PathIndex last_index,
                          const uint8_t *query, size_t path_nbytes);

// Add a PackedPath, using a FileFilter to reduce to a subset of colours
// `find` Specifies if we should try to find a duplicate first
// Returns PATH_NULL if no colours set in colour subset
PathIndex path_store_find_or_add(PathStore *store,
                                 BinaryKmer bkmer, PathIndex last_index,
                                 const uint8_t *packed, size_t path_nbytes,
                                 const FileFilter *fltr, bool find,
                                 bool *added);

// If compatible, a FileFilter can be read straight into a PathStore without
// parsing each path, one-by-one (much faster!)
#define path_store_fltr_compatible(st,fltr) \
        ((fltr)->nofilter && \
         roundup_bits2bytes((fltr)->filencols) == (st)->colset_bytes)

//
// Print
//
void path_store_print(const PathStore *pstore);

void path_store_print_path(const PathStore *paths, PathIndex index);
void path_store_print_all(const PathStore *paths);

// packed points to <PathLen><PackedSeq>
void print_path(hkey_t hkey, const uint8_t *packed, const PathStore *pstore);

//
// Data checks for debugging / testing
//
bool path_store_data_integrity_check(const uint8_t *data, size_t size,
                                        size_t colbytes);

bool path_store_integrity_check(const PathStore *pstore);

#endif /* BINARY_PATH_H_ */
