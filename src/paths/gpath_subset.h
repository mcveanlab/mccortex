#ifndef GPATH_SUBSET_H_
#define GPATH_SUBSET_H_

#include "gpath.h"
#include "gpath_set.h"

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(gpath_ptr_buf, GPathPtrBuffer, GPath*);

typedef struct
{
  GPathSet *gpset;
  GPathPtrBuffer list; // Sort and remove entries from this array
  bool is_sorted;
} GPathSubset;

void gpath_subset_alloc(GPathSubset *subset);
void gpath_subset_dealloc(GPathSubset *subset);
void gpath_subset_reset(GPathSubset *subset);
void gpath_subset_init(GPathSubset *subset, GPathSet *gpset);
void gpath_subset_add(GPathSubset *subset, GPath *path);
void gpath_subset_sort(GPathSubset *subset);

// Load linked list pointed to by first. If NULL load nothing.
// Does not reset subset before loading
void gpath_subset_load_llist(GPathSubset *subset, GPath *first);

// Load all paths from a given set
// Does not reset subset before loading
void gpath_subset_load_set(GPathSubset *subset);

// Update the linked list of paths in set `subset->set`
// You still need to update the pointer to the first item (e.g. in GPathStore)
void gpath_subset_update_linkedlist(GPathSubset *subset);

/**
 * Remove duplicate entries e.g.
 *  {T,TT,TT} -> {T,TT}
 */
void gpath_subset_rmdup(GPathSubset *subset);

/**
 * Remove redundant entries such as duplicates and substrings e.g.
 *  {T,TT,TT} -> {TT}
 *  {A,C,CG,CGC} -> {A,CGC}
 */
void gpath_subset_rmsubstr(GPathSubset *subset);

/**
 * Remove entries from `src` that are in `dst`, copying over sample counts
 */
void gpath_subset_merge(GPathSubset *dst, GPathSubset *src);

#endif /* GPATH_SUBSET_H_ */
