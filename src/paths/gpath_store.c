#include "global.h"
#include "gpath_store.h"


// @split_linked_lists whether you intend to have traverse linked list and
//                     all linked list separate
size_t gpath_store_min_mem(size_t graph_capacity, bool split_linked_lists)
{
  return graph_capacity*sizeof(GPath*) * (split_linked_lists ? 2 : 1);
}

void gpath_store_alloc(GPathStore *gpstore, size_t ncols, size_t graph_capacity,
                       size_t mem, bool count_nseen, bool split_linked_lists)
{
  memset(gpstore, 0, sizeof(*gpstore));

  size_t min_mem = gpath_store_min_mem(count_nseen, split_linked_lists);
  ctx_assert2(min_mem < mem, "%zu is not less than %zu", min_mem, mem);

  // Don't allow resizing - this allows use to the thread safe
  // and control memory usage (die when we fill up allowed memory)
  gpath_set_alloc(&gpstore->gpset, ncols, mem - min_mem, false, count_nseen);

  gpstore->graph_capacity = graph_capacity;
  gpstore->paths_all = ctx_calloc(graph_capacity, sizeof(GPath*));
  gpstore->paths_traverse = gpstore->paths_all;
  // gpstore->kmer_locks = ctx_calloc((graph_capacity+7)/8, 1);
}

void gpath_store_dealloc(GPathStore *gpstore)
{
  gpath_set_dealloc(&gpstore->gpset);
  gpath_store_merge_read_write(gpstore);
  ctx_free(gpstore->paths_all);
  // ctx_free(gpstore->kmer_locks);
  memset(gpstore, 0, sizeof(*gpstore));
}

void gpath_store_reset(GPathStore *gpstore)
{
  gpath_set_reset(&gpstore->gpset);
  gpstore->num_kmers_with_paths = gpstore->num_paths = gpstore->path_bytes = 0;
  memset(gpstore->paths_all, 0, gpstore->graph_capacity * sizeof(GPath*));
  if(gpstore->paths_traverse != gpstore->paths_all)
    memset(gpstore->paths_traverse, 0, gpstore->graph_capacity * sizeof(GPath*));
}

void gpath_store_split_read_write(GPathStore *gpstore)
{
  if(gpstore->num_paths > 0) {
    // Copy current paths over to path set to be updated
    status("[PathStore] Creating separate read/write GraphPath linked lists");
    size_t mem = gpstore->graph_capacity * sizeof(GPath*);
    gpstore->paths_traverse = ctx_calloc(gpstore->graph_capacity, sizeof(GPath*));
    memcpy(gpstore->paths_traverse, gpstore->paths_all, mem);
  }
}

void gpath_store_merge_read_write(GPathStore *gpstore)
{
  if(gpstore->paths_traverse != gpstore->paths_all)
  {
    status("[PathStore] Merging read/write GraphPath linked lists");
    ctx_free(gpstore->paths_traverse);
    gpstore->paths_traverse = gpstore->paths_all;
  }
}

GPath* gpath_store_fetch(const GPathStore *gpstore, hkey_t hkey)
{
  return gpstore->paths_all[hkey];
}

GPath* gpath_store_fetch_traverse(const GPathStore *gpstore, hkey_t hkey)
{
  return gpstore->paths_traverse[hkey];
}

// Update stats after removing a path
void gpstore_path_removal_update_stats(GPathStore *gpstore, GPath *gpath)
{
  // Update stats
  size_t nbytes = (gpath->num_juncs+3)/4;
  __sync_fetch_and_sub((volatile uint64_t*)&gpstore->num_paths, 1);
  __sync_fetch_and_sub((volatile uint64_t*)&gpstore->path_bytes, nbytes);
}

// You do not need to acquire the kmer lock before calling this function
static void _gpstore_add_to_llist_mt(GPathStore *gpstore, hkey_t hkey, GPath *gpath)
{
  // Update stats
  size_t nbytes = (gpath->num_juncs+3)/4;
  size_t new_kmer = (gpstore->paths_all[hkey] == NULL);
  __sync_fetch_and_add((volatile uint64_t*)&gpstore->num_kmers_with_paths, new_kmer);
  __sync_fetch_and_add((volatile uint64_t*)&gpstore->num_paths, 1);
  __sync_fetch_and_add((volatile uint64_t*)&gpstore->path_bytes, nbytes);

  // Add to linked list
  do {
    gpath->next = *(GPath *volatile const*)&gpstore->paths_all[hkey];
  }
  while(!__sync_bool_compare_and_swap((GPath *volatile*)&gpstore->paths_all[hkey],
                                      gpath->next, gpath));
}

// Linear search to find a given path
GPath* gpstore_find(const GPathStore *gpstore, hkey_t hkey, GPathNew find)
{
  GPath *gpath = gpath_store_fetch(gpstore, hkey);
  for(; gpath != NULL; gpath = gpath->next)
    if(gpaths_are_equal(*gpath, find))
      return gpath;
  return NULL;
}

// Always adds
// Note: it is not safe to call _add and _find_add simultaneously, since _add
//       avoids the use of locks.
GPath* gpath_store_add_mt(GPathStore *gpstore, hkey_t hkey, GPathNew newgpath)
{
  ctx_assert(newgpath.seq != NULL);

  GPath *gpath = gpath_set_add_mt(&gpstore->gpset, newgpath);
  _gpstore_add_to_llist_mt(gpstore, hkey, gpath);

  return gpath;
}


// NOT USED ATM
// Linear time search to find or add a given path
// Note: it is not safe to call _add and _find_add simultaneously, since _add
//       avoids the use of locks.
GPath* gpath_store_find_add_mt(GPathStore *gpstore,
                               hkey_t hkey, GPathNew newgpath,
                               bool *found)
{
  ctx_assert(newgpath.seq != NULL);

  GPath *gpath;
  *found = true;

  // Get lock for kmer
  bitlock_yield_acquire(gpstore->kmer_locks, hkey);

  // Add if not found
  if((gpath = gpstore_find(gpstore, hkey, newgpath)) == NULL) {
    gpath = gpath_set_add_mt(&gpstore->gpset, newgpath);
    *found = false;
  }

  // Release kmer lock
  bitlock_release(gpstore->kmer_locks, hkey);

  _gpstore_add_to_llist_mt(gpstore, hkey, gpath);

  return gpath;
}
