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

// Remove duplicate entries
// {T,TT,TT} -> {T,TT}
void gpath_subset_rmdup(GPathSubset *subset)
{
  if(subset->list.len <= 1) return;
  if(!subset->is_sorted) gpath_subset_sort(subset);

  size_t i, j, len = subset->list.len, ncols = subset->gpset->ncols;
  GPath **list = subset->list.b;

  for(i = 0, j = 1; j < len; j++) {
    if(binary_seqs_cmp(list[i]->seq, list[i]->num_juncs,
                       list[j]->seq, list[j]->num_juncs) == 0)
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

// Remove redundant entries such as duplicates and substrings
// {T,TT,TT} -> {TT}
// {A,C,CG,CGC} -> {A,CGC}
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

          // j can't be a subset of i if it's longer
          if(list[j]->num_juncs > list[i]->num_juncs ||
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
  for(i = j = 0; i < len; i++) {
    if(list[i] != NULL)
      list[j++] = list[i];
  }
  subset->list.len = j;
}

// Remove entries from `src` that are in `dst`,
// if `rmsubstr` then also remove substring matches
// Note: does not remove duplicates from `dst`,
//       call gpath_subset_rmsubstr() to do that
// We add colours to `dst` if they are removed from `src`
//
// Thread safe for accesses to `dst` but not `src`. In other words, this thread
// must be the only one reading/writing `src`, but multiple threads can be
// accessing `dst`
void gpath_subset_merge(GPathSubset *dst, GPathSubset *src, bool rmsubstr)
{
  ctx_assert2(dst->gpset->ncols == src->gpset->ncols, "%zu vs %zu",
              dst->gpset->ncols, src->gpset->ncols);

  gpath_ptr_buf_capacity(&dst->list, dst->list.len+src->list.len);
  size_t dstlen = dst->list.len, srclen = src->list.len;
  size_t i = 0, j = 0, k, ncols = dst->gpset->ncols, min_num_juncs;
  int cmp;

  GPath **dstlist = dst->list.b;
  GPath **srclist = src->list.b;

  if(dstlen == 0 || srclen == 0) return;

  while(i < dstlen && j < srclen)
  {
    if(dstlist[i]->num_juncs < srclist[j]->num_juncs) {
      // No way that srclist[j] can be a substr of dstlist[i]
      cmp = binary_seqs_cmp(dstlist[i]->seq, dstlist[i]->num_juncs,
                            srclist[j]->seq, srclist[j]->num_juncs);
    }
    else if(!rmsubstr) {
      cmp = binary_seqs_cmp(dstlist[i]->seq, dstlist[i]->num_juncs,
                            srclist[j]->seq, srclist[j]->num_juncs);
      if(dstlist[i]->num_juncs == srclist[j]->num_juncs && cmp == 0) {
        // paths match, steal colours and remove it
        gpath_colset_or_mt(dstlist[i], srclist[j], ncols);
        gpath_set_nseen_sum_mt(dstlist[i], dst->gpset,
                               srclist[j], src->gpset);
        srclist[j] = NULL;
      }
    }
    else {
      // Find paths that srclist[j] is a subset of
      // Remove those colours from srclist[j]
      // Keep checking until colset is zero
      for(k = i; k < dstlen; k++)
      {
        min_num_juncs = MIN2(dstlist[k]->num_juncs, srclist[j]->num_juncs);

        cmp = binary_seqs_cmp(dstlist[i]->seq, min_num_juncs,
                              srclist[j]->seq, min_num_juncs);
        if(cmp != 0) break;

        if(dstlist[k]->num_juncs == srclist[j]->num_juncs) {
          // paths match, steal colours and remove it
          gpath_colset_or_mt(dstlist[i], srclist[j], ncols);
          gpath_set_nseen_sum_mt(dstlist[i], dst->gpset,
                                 srclist[j], src->gpset);
          srclist[j] = NULL;
          break;
        }

        // srclist[j] is substring of dstlist[k]
        // remove colours from j that are in k
        if(gpath_colset_rm_intersect(dstlist[k], srclist[j], ncols) == 0) {
          srclist[j] = NULL;
          break;
        }
      }
    }

    if(cmp < 0) i++;
    else j++;
  }

  // Remove NULLs from src
  for(i = j = 0; i < src->list.len; i++)
    if(srclist[i] != NULL)
      srclist[j++] = srclist[i];

  src->list.len = j;
}

