#include "global.h"
#include "sorted_path_set.h"
#include "binary_seq.h"

//
// SortedPathSet allows fast comparison of sets of paths and duplicate removal
//

void sorted_path_set_alloc(SortedPathSet *set)
{
  bytebuf_alloc(&set->seqs, 512);
  pentrybuf_alloc(&set->members, 512);
}

void sorted_path_set_dealloc(SortedPathSet *set)
{
  bytebuf_dealloc(&set->seqs);
  pentrybuf_dealloc(&set->members);
}

// Sort in orientation (FORWARD, then REVERSE)
// Then sort in lexicographic order
//   a
//   ab
//   ac
//   b
//   bc
//   bcd
static int _sorted_path_entry_cmp(const void *aa, const void *bb)
{
  const PathEntry *a = (const PathEntry*)aa, *b = (const PathEntry*)bb;
  int ret;
  if((ret = (int)a->orient - b->orient) != 0) return ret;
  if((ret = binary_seqs_cmp(a->seq, a->plen, b->seq, b->plen)) != 0) return ret;
  return a->pindex - b->pindex;
}

// Load a PathSet from a given PathStore
// SortedPathSet allows fast comparison of sets and duplicate removal
//  `cbytes` is number of bytes in colourset of both PathStore and SortedPathSet
void sorted_path_set_init2(SortedPathSet *set, hkey_t hkey, size_t cbytes,
                           const uint8_t *store, PathIndex pindex,
                           const FileFilter *fltr)
{
  set->hkey = hkey;
  set->cbytes = cbytes;

  ByteBuffer *bytebuf = &set->seqs;
  PathEntryBuffer *membuf = &set->members;

  bytebuf_reset(bytebuf);
  pentrybuf_reset(membuf);

  // Fetch paths into sorted set
  PathLen plen = 0;
  Orientation orient = FORWARD;
  size_t i, nbytes, intocol, fromcol;

  while(pindex != PATH_NULL)
  {
    const uint8_t *packed = store+pindex;
    const uint8_t *store_colset = packedpath_get_colset(packed);
    const uint8_t *store_seq = packedpath_seq(packed, cbytes);

    packedpath_get_len_orient(packed, cbytes, &plen, &orient);

    // Add entry
    PathEntry entry = {.seq = NULL, .pindex = pindex,
                       .orient = orient, .plen = plen};
    pentrybuf_add(membuf, entry);

    // Copy colset+sequence to buffer
    nbytes = (plen+3)/4;
    bytebuf_ensure_capacity(bytebuf, bytebuf->len + cbytes + nbytes);

    uint8_t *set_colset = bytebuf->data + bytebuf->len;
    uint8_t *set_seq = set_colset + cbytes;

    // Copy colour bitset into sorted set
    if(fltr == NULL) {
      memcpy(set_colset, store_colset, cbytes);
    } else {
      // Clear memory for colour bitset
      memset(set_colset, 0, cbytes);

      // Copy over bitset, one bit at a time
      for(i = 0; i < fltr->ncols; i++) {
        intocol = file_filter_intocol(fltr, i);
        fromcol = file_filter_fromcol(fltr, i);
        bitset_cpy(set_colset, intocol, bitset_get(store_colset, fromcol));
      }
    }

    // Copy sequence into sorted set
    memcpy(set_seq, store_seq, nbytes);
    bytebuf->len += cbytes + nbytes;

    pindex = packedpath_get_prev(packed);
  }

  // Set pointers for paths
  uint8_t *seq = set->seqs.data;

  for(i = 0; i < set->members.len; i++)
  {
    membuf->data[i].seq = seq + cbytes;
    seq += cbytes + (membuf->data[i].plen + 3)/4;
  }

  // Sort paths
  qsort(set->members.data, set->members.len, sizeof(PathEntry),
        _sorted_path_entry_cmp);
}

void sorted_path_set_init(SortedPathSet *set, const PathStore *ps, hkey_t hkey)
{
  PathIndex pindex = pstore_get_pindex(ps, hkey);
  sorted_path_set_init2(set, hkey, ps->colset_bytes, ps->store, pindex, NULL);
}

// Remove redundant entries such as duplicates and substrings
// {T,TT,TT} -> {TT}
// {A,C,CG,CGC} -> {A,CGC}
void sorted_path_set_slim(SortedPathSet *set)
{
  const size_t num_members = set->members.len;
  if(num_members == 0) return;

  size_t i, j;
  bool potential_empty_sets = false;
  PathEntry *members = set->members.data;
  PathLen min_plen;
  uint8_t *cset_i, *cset_j; // colour sets

  // Work backwards to remove colours from subsumed paths
  for(i = num_members-1; i > 0; i--)
  {
    cset_i = sorted_path_colset(&members[i],set);

    // If sample has no colours skip it
    if(!packedpath_is_colset_zero(cset_i, set->cbytes))
    {
      // work backwards over paths as j as they match
      for(j = i-1; j != SIZE_MAX && members[j].plen <= members[i].plen; j--)
      {
        min_plen = MIN2(members[i].plen, members[j].plen);
        if(binary_seqs_cmp(members[i].seq, min_plen, members[j].seq, min_plen) != 0) {
          break;
        }

        cset_j = sorted_path_colset(&members[j],set);

        if(members[i].plen == members[j].plen) {
          // paths i,j match, steal colours, zero it
          packedpath_cpy_zero_colsets(cset_i, cset_j, set->cbytes);
        }
        else {
          // path j is substring of i, remove colours from j that are in i
          packedpath_colset_rm_intersect(cset_j, cset_i, set->cbytes);
        }

        potential_empty_sets = true;
      }
    }
  }

  if(!potential_empty_sets) return;

  // loop over entries and remove empty ones
  for(i = j = 0; i < num_members; i++) {
    cset_i = sorted_path_colset(&members[i],set);
    if(!packedpath_is_colset_zero(cset_i, set->cbytes)) {
      members[j++] = members[i];
    }
  }

  set->members.len = j;
}

// Updates entry0
static inline void _update_store(const PathEntry *entry0, const SortedPathSet *set0,
                                 const PathEntry *entry1, const SortedPathSet *set1,
                                 uint8_t *store)
{
  if(store == NULL) return;
  ctx_assert(set0->cbytes == set1->cbytes);

  uint8_t *cset0 = sorted_path_colset(entry0, set0);
  const uint8_t *cset1 = sorted_path_colset(entry1, set1);

  // packedpath_colsets_or returns true if colset was changed
  if(packedpath_colsets_or(cset0, cset1, set0->cbytes))
  {
    // copy changes to store
    uint8_t *ptr = packedpath_get_colset(store + entry0->pindex);
    memcpy(ptr, cset0, set0->cbytes);
  }
}

// Remove entries from set that are in the filter set
// if `rmsubstr` then also remove substring matches
// Note: does not remove duplicates from `set`,
//       call sorted_path_set_slim() to do that
// if passed store, we use it to update colourset of paths in out_set
void sorted_path_set_merge(SortedPathSet *out_set, SortedPathSet *in_set,
                           bool rmsubstr, uint8_t *store)
{
  // for each entry in in_set, run through out_set and check if its contained

  const size_t num_out = out_set->members.len, num_in = in_set->members.len;
  const size_t cbytes = in_set->cbytes;
  size_t i = 0, j = 0, k = 0, x = 0;
  PathLen min_plen;
  int cmp;

  ctx_assert(in_set->cbytes == out_set->cbytes);
  ctx_assert(in_set != out_set);

  if(num_in == 0 || num_out == 0) return;

  PathEntry *list0 = out_set->members.data;
  PathEntry *list1 = in_set->members.data;
  uint8_t *cset_in, *cset_out; // coloursets

  if(rmsubstr)
  {
    while(i < num_out && j < num_in) {
      cset_in = sorted_path_colset(&list1[j], in_set);

      if(list0[i].plen < list1[j].plen) {
        // No way that list1[j] can be a substr of list0[i]
        cmp = binary_seqs_cmp(list0[i].seq, list0[i].plen,
                              list1[j].seq, list1[j].plen);
      }
      else {
        // Find paths that list1[j] is a subset of
        // Remove those colours from list1[j]
        // Keep checking until colset is zero
        for(k = i; k < num_out; k++)
        {
          min_plen = MIN2(list0[k].plen, list1[j].plen);
          cmp = binary_seqs_cmp(list0[k].seq, min_plen, list1[j].seq, min_plen);
          if(cmp != 0) break;

          // list1[j] is substring of or match to list0[k]

          if(list0[k].plen == list1[j].plen) {
            // paths match - add colour membership to list0[k]
            _update_store(&list0[k], out_set, &list1[j], in_set, store);
          }

          // remove colours from j that are in k
          cset_out = sorted_path_colset(&list0[k], in_set);
          packedpath_colset_rm_intersect(cset_in, cset_out, cbytes);
          if(packedpath_is_colset_zero(cset_in, cbytes)) break;
        }
      }

      if(cmp < 0) i++;
      else if(cmp > 0 && !packedpath_is_colset_zero(cset_in, cbytes)) {
        list0[x++] = list0[j++];
      } else j++;
    }
  }
  else
  {
    while(i < num_out && j < num_in) {
      cmp = binary_seqs_cmp(list0[i].seq, list0[i].plen, list1[j].seq, list1[j].plen);
      if(cmp < 0) i++;
      else if(cmp == 0) {
        _update_store(&list0[i], out_set, &list1[j], in_set, store); j++;
      } else list0[x++] = list0[j++];
    }
  }

  for(; j < num_in; j++) {
    cset_in = sorted_path_colset(&list1[j], in_set);
    if(!packedpath_is_colset_zero(cset_in, cbytes)) list1[x++] = list1[j];
  }

  in_set->members.len = x;
}
