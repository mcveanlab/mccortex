#include "global.h"
#include "gpath_subset.h"
#include "binary_seq.h"

void gpath_subset_alloc(GPathSubset *subset)
{
  gpath_ptr_buf_alloc(&subset->list, 16);
  gpath_subset_reset(subset);
}

void gpath_subset_dealloc(GPathSubset *subset)
{
  gpath_ptr_buf_dealloc(&subset->list);
  memset(subset, 0, sizeof(GPathSubset));
}

void gpath_subset_reset(GPathSubset *subset)
{
  gpath_ptr_buf_reset(&subset->list);
  subset->is_sorted = false;
}

void gpath_subset_init(GPathSubset *subset, GPathSet *gpset)
{
  gpath_subset_reset(subset);
  subset->gpset = gpset;
}

void gpath_subset_add(GPathSubset *subset, GPath *path)
{
  gpath_ptr_buf_add(&subset->list, path);
  subset->is_sorted = false;
}

void gpath_subset_sort(GPathSubset *subset)
{
  qsort(subset->list.b, subset->list.len, sizeof(GPath*), gpath_cmp_void);
  subset->is_sorted = true;
}

// Load linked list pointed to by first. If NULL load nothing.
// Does not reset subset before loading
void gpath_subset_load_llist(GPathSubset *subset, GPath *first)
{
  while(first != NULL)
  {
    gpath_ptr_buf_add(&subset->list, first);
    first = first->next;
  }
}

// Load all paths from a given set
// Does not reset subset before loading
void gpath_subset_load_set(GPathSubset *subset)
{
  GPathSet *gpset = subset->gpset;
  size_t i;
  for(i = 0; i < gpset->entries.len; i++)
    gpath_ptr_buf_add(&subset->list, &gpset->entries.b[i]);
}

// Update the linked list of paths in set `subset->gpset`
// You still need to update the pointer to the first item in GPathStore
void gpath_subset_update_linkedlist(GPathSubset *subset)
{
  if(subset->list.len == 0) return;
  size_t i;
  for(i = 0; i+1 < subset->list.len; i++)
    subset->list.b[i]->next = subset->list.b[i+1];
  subset->list.b[subset->list.len-1]->next = NULL;
}

/**
 * Remove duplicate entries e.g.
 *  {T,TT,TT} -> {T,TT}
 */
void gpath_subset_rmdup(GPathSubset *subset)
{
  if(subset->list.len <= 1) return;
  if(!subset->is_sorted) gpath_subset_sort(subset);

  size_t i, j, len = subset->list.len, ncols = subset->gpset->ncols;
  GPath **list = subset->list.b;

  for(i = 0, j = 1; j < len; j++) {
    if(gpath_cmp(list[i], list[j]) == 0)
    {
      gpath_colset_or_mt(list[i], list[j], ncols);
      gpath_set_nseen_sum_mt(list[i], subset->gpset,
                             list[j], subset->gpset);
    }
    else {
      list[++i] = list[j];
    }
  }

  // i is the index of the last unique entry
  subset->list.len = i+1;
}

/**
 * Remove redundant entries such as duplicates and substrings e.g.
 *  {T,TT,TT} -> {TT}
 *  {A,C,CG,CGC} -> {A,CGC}
 */
void gpath_subset_rmsubstr(GPathSubset *subset)
{
  if(subset->list.len <= 1) return;
  if(!subset->is_sorted) gpath_subset_sort(subset);

  size_t i, j, len = subset->list.len, min_juncs, ncols = subset->gpset->ncols;
  GPath **list = subset->list.b;

  // Work backwards to remove colours from subsumed paths
  for(i = len-1; i > 0; i--) {
    if(list[i] != NULL) {
      // work backwards over paths as j as they match
      for(j = i-1; j != SIZE_MAX; j--)
      {
        if(list[j] != NULL)
        {
          min_juncs = MIN2(list[i]->num_juncs, list[j]->num_juncs);

          // j can't be a subset of i if it's longer,
          // or orientations don't match
          if(list[i]->num_juncs < list[j]->num_juncs ||
             list[i]->orient != list[j]->orient ||
             binary_seqs_cmp(list[i]->seq, min_juncs,
                             list[j]->seq, min_juncs) != 0)
          {
            break;
          }
          else if(list[i]->num_juncs == list[j]->num_juncs)
          {
            // paths match, steal colours from j and remove it
            gpath_colset_or_mt(list[i], list[j], ncols);
            gpath_set_nseen_sum_mt(list[i], subset->gpset,
                                   list[j], subset->gpset);
            list[j] = NULL;
          }
          else
          {
            // path j is substring of i, remove colours from j that are in i
            // then remove path j only if all colours removed
            if(gpath_colset_rm_intersect(list[i], list[j], ncols) == 0)
              list[j] = NULL;
          }
        }
      }
    }
  }

  // loop over entries and remove empty ones
  for(i = j = 0; i < len; i++)
    if(list[i] != NULL)
      list[j++] = list[i];

  subset->list.len = j;
}

/**
 * Remove entries from `src` that are in `dst`, copying over sample counts
 */
void gpath_subset_merge(GPathSubset *dst, GPathSubset *src)
{
  ctx_assert2(dst->gpset->ncols == src->gpset->ncols, "%zu vs %zu",
              dst->gpset->ncols, src->gpset->ncols);

  if(!dst->is_sorted) gpath_subset_sort(dst);
  if(!src->is_sorted) gpath_subset_sort(src);

  size_t i = 0, j = 0, ncols = dst->gpset->ncols;
  int cmp;

  GPath **dstlist = dst->list.b;
  GPath **srclist = src->list.b;

  if(dst->list.len == 0 || src->list.len == 0) return;

  while(i < dst->list.len && j < src->list.len)
  {
    cmp = gpath_cmp(dstlist[i], srclist[j]);

    if(cmp < 0) i++;
    else if(cmp > 0) j++;
    else {
      // paths match, steal colours and remove it
      gpath_colset_or_mt(dstlist[i], srclist[j], ncols);
      gpath_set_nseen_sum_mt(dstlist[i], dst->gpset,
                             srclist[j], src->gpset);
      srclist[j] = NULL;
      j++;
    }
  }

  // Remove NULLs from src
  for(i = j = 0; i < src->list.len; i++)
    if(srclist[i] != NULL)
      srclist[j++] = srclist[i];

  src->list.len = j;
}

