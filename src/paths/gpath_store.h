#ifndef GPATH_STORE_H_
#define GPATH_STORE_H_

#include "gpath_set.h"

// GPathStore is a map from {[kmer/hkey] -> [GPath linked list]}
typedef struct
{
  uint64_t num_kmers_with_paths, num_paths, path_bytes;
  uint64_t graph_capacity;
  GPathSet gpset;
  GPath **paths_all, **paths_traverse;
} GPathStore;

// @split_linked_lists whether you intend to have traverse linked list and
//                     all linked list separate
//                     (i.e call gpath_store_split_read_write())
size_t gpath_store_mem(size_t graph_capacity, bool split_linked_lists);

// If num_paths != 0, we ensure at least num_paths capacity
void gpath_store_alloc(GPathStore *gpstore, size_t ncols, size_t graph_capacity,
                       size_t num_paths, size_t mem,
                       bool count_nseen, bool split_linked_lists);

void gpath_store_dealloc(GPathStore *gpstore);
void gpath_store_reset(GPathStore *gpstore);

void gpath_store_split_read_write(GPathStore *gpstore);
void gpath_store_merge_read_write(GPathStore *gpstore);

// Traversal paths are a subset of all paths
#define gpath_store_use_traverse(gpstore) ((gpstore)->paths_traverse != NULL)
GPath* gpath_store_fetch(const GPathStore *gpstore, hkey_t hkey);
GPath* gpath_store_fetch_traverse(const GPathStore *gpstore, hkey_t hkey);

GPath* gpstore_find(const GPathStore *gpstore, hkey_t hkey, GPathNew find);

// Always adds
// colset sequence will be zeroed on return
GPath* gpath_store_add_mt(GPathStore *gpstore, hkey_t hkey, GPathNew newgpath);

/*
// Not used
// Linear time search to find or add a given path
// if added, colset will be zero'd and gpath.seq will be copied
GPath* gpath_store_find_add_mt(GPathStore *gpstore,
                               hkey_t hkey, GPathNew newgpath,
                               bool *found);
*/

// Update stats after removing a path
void gpstore_path_removal_update_stats(GPathStore *gpstore, GPath *gpath);

#endif /* GPATH_STORE_H_ */
