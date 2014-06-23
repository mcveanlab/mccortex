#include "global.h"
#include "path_set.h"
#include "binary_seq.h"
#include "sort_r/sort_r.h"

//
// PathSet allows fast comparison of sets of paths and duplicate removal
//

void path_set_alloc(PathSet *set)
{
  byte_buf_alloc(&set->seqs, 512);
  pentrybuf_alloc(&set->members, 512);
  path_set_reset(set);
}

void path_set_dealloc(PathSet *set)
{
  byte_buf_dealloc(&set->seqs);
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
static int _path_entry_cmp(const void *aa, const void *bb, void *ptr)
{
  const PathEntry *a = (const PathEntry*)aa, *b = (const PathEntry*)bb;
  const PathSet *set = (const PathSet*)ptr;
  int ret;
  if((ret = (int)a->orient - b->orient) != 0) return ret;
  ret = binary_seqs_cmp(path_set_seq(a,set), a->plen,
                        path_set_seq(b,set), b->plen);
  if(ret != 0) return ret;
  return a->pindex - b->pindex;
}

void path_set_sort(PathSet *set)
{
  // Sort paths
  sort_r(set->members.data, set->members.len, sizeof(PathEntry),
         _path_entry_cmp, set);
  set->sorted = true;
}

// Empty the path set
void path_set_reset(PathSet *set)
{
  byte_buf_reset(&set->seqs);
  pentrybuf_reset(&set->members);
  set->sorted = false;
  set->cbytes = 0;
}

// Load a PathSet from a given PathStore
// PathSet allows fast comparison of sets and duplicate removal
//  `cbytes` is number of bytes in colourset of both PathStore and PathSet
void path_set_init2(PathSet *set, size_t cbytes,
                    const uint8_t *store, PathIndex pindex,
                    const FileFilter *fltr)
{
  path_set_reset(set);
  path_set_load(set, cbytes, store, pindex, fltr);
}

void path_set_init(PathSet *set, const PathStore *ps, hkey_t hkey)
{
  PathIndex pindex = pstore_get_pindex(ps, hkey);
  path_set_init2(set, ps->colset_bytes, ps->store, pindex, NULL);
}

// Load paths into the set. If not reset, appends.
void path_set_load(PathSet *set, size_t cbytes,
                   const uint8_t *store, PathIndex pindex,
                   const FileFilter *fltr)
{
  ctx_assert(set->cbytes == 0 || cbytes == set->cbytes);

  set->cbytes = cbytes;

  ByteBuffer *bytebuf = &set->seqs;
  PathEntryBuffer *membuf = &set->members;

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
    PathEntry entry = {.seq = bytebuf->len + cbytes, .pindex = pindex,
                       .orient = orient, .plen = plen, .count = 0};
    pentrybuf_add(membuf, entry);

    // Copy colset+sequence to buffer
    nbytes = (plen+3)/4;
    byte_buf_ensure_capacity(bytebuf, bytebuf->len + cbytes + nbytes);

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

  set->sorted = false;
}

// Remove redundant entries such as duplicates, substrings and
// paths with no colours
// {T,TT,TT} -> {TT}
// {A,C,CG,CGC} -> {A,CGC}
void path_set_slim(PathSet *set)
{
  const size_t num_members = set->members.len;
  if(num_members == 0) return;

  if(!set->sorted) path_set_sort(set);

  size_t i, j;
  bool potential_empty_sets = false;
  PathEntry *members = set->members.data;
  PathLen min_plen;
  uint8_t *cset_i, *cset_j; // colour sets
  const uint8_t *seq_i, *seq_j;

  // Work backwards to remove colours from subsumed paths
  for(i = num_members-1; i > 0; i--)
  {
    seq_i = path_set_seq(&members[i],set);
    cset_i = path_set_colset(&members[i],set);

    // If sample has no colours skip it
    if(packedpath_is_colset_zero(cset_i, set->cbytes)) {
      potential_empty_sets = true;
    }
    else {
      // work backwards over paths as j as they match
      for(j = i-1; j != SIZE_MAX && members[j].plen <= members[i].plen; j--)
      {
        min_plen = MIN2(members[i].plen, members[j].plen);

        seq_j = path_set_seq(&members[j],set);
        cset_j = path_set_colset(&members[j],set);

        if(binary_seqs_cmp(seq_i, min_plen, seq_j, min_plen) != 0)
          break;

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

  // We haven't checked the first path entry
  if(!potential_empty_sets) {
    cset_i = path_set_colset(&members[0],set);
    potential_empty_sets = packedpath_is_colset_zero(cset_i, set->cbytes);
  }

  if(!potential_empty_sets) return;

  // loop over entries and remove empty ones
  for(i = j = 0; i < num_members; i++) {
    cset_i = path_set_colset(&members[i],set);
    if(!packedpath_is_colset_zero(cset_i, set->cbytes)) {
      members[j++] = members[i];
    }
  }

  set->members.len = j;
}

// Updates entry0
static inline void _update_store(const PathEntry *entry0, const PathSet *set0,
                                 const PathEntry *entry1, const PathSet *set1,
                                 uint8_t *store)
{
  if(store == NULL) return;
  ctx_assert(set0->cbytes == set1->cbytes);

  uint8_t *cset0 = path_set_colset(entry0, set0);
  const uint8_t *cset1 = path_set_colset(entry1, set1);

  // packedpath_colsets_or returns true if colset was changed
  if(packedpath_colsets_or(cset0, cset1, set0->cbytes))
  {
    // copy changes to store
    uint8_t *ptr = packedpath_get_colset(store + entry0->pindex);
    memcpy(ptr, cset0, set0->cbytes);
  }
}

// DEV: does not update number of path bytes - doesn't need to!
// Remove entries from set that are in the filter set
// if `rmsubstr` then also remove substring matches
// Note: does not remove duplicates from `out_set`,
//       call path_set_slim() to do that
// if passed store, we use it to update colourset of paths in out_set
void path_set_merge(PathSet *out_set, PathSet *in_set,
                    bool rmsubstr, uint8_t *store)
{
  // for each entry in in_set, run through out_set and check if its contained

  if(!out_set->sorted) path_set_sort(out_set);
  if(!in_set->sorted) path_set_sort(in_set);

  const size_t num_out = out_set->members.len, num_in = in_set->members.len;
  const size_t cbytes = in_set->cbytes;
  size_t i = 0, j = 0, k = 0, x = 0;
  PathLen min_plen;
  int cmp;

  ctx_assert(in_set->cbytes == out_set->cbytes);
  ctx_assert(in_set != out_set);

  if(num_in == 0 || num_out == 0) return;

  PathEntry *outlist = out_set->members.data;
  PathEntry *inlist = in_set->members.data;
  uint8_t *cset_in, *cset_out; // coloursets
  const uint8_t *outseq, *inseq, *kseq;

  if(rmsubstr)
  {
    while(i < num_out && j < num_in) {
      cset_in = path_set_colset(&inlist[j], in_set);
      outseq = path_set_seq(&outlist[i],out_set);
      inseq = path_set_seq(&inlist[j],in_set);

      if(outlist[i].plen < inlist[j].plen) {
        // No way that inlist[j] can be a substr of outlist[i]
        cmp = binary_seqs_cmp(outseq, outlist[i].plen, inseq, inlist[j].plen);
      }
      else {
        // Find paths that inlist[j] is a subset of
        // Remove those colours from inlist[j]
        // Keep checking until colset is zero
        for(k = i; k < num_out; k++)
        {
          kseq = path_set_seq(&outlist[k],out_set);
          min_plen = MIN2(outlist[k].plen, inlist[j].plen);
          cmp = binary_seqs_cmp(kseq, min_plen, inseq, min_plen);
          if(cmp != 0) break;

          // inlist[j] is substring of or match to outlist[k]

          if(outlist[k].plen == inlist[j].plen) {
            // paths match - add colour membership to outlist[k]
            _update_store(&outlist[k], out_set, &inlist[j], in_set, store);
          }

          // remove colours from j that are in k
          cset_out = path_set_colset(&outlist[k], in_set);
          packedpath_colset_rm_intersect(cset_in, cset_out, cbytes);
          if(packedpath_is_colset_zero(cset_in, cbytes)) break;
        }
      }

      if(cmp < 0) i++;
      else if(cmp > 0 && !packedpath_is_colset_zero(cset_in, cbytes)) {
        outlist[x++] = outlist[j++];
      } else j++;
    }
  }
  else
  {
    while(i < num_out && j < num_in) {
      outseq = path_set_seq(&outlist[i],out_set);
      inseq = path_set_seq(&inlist[j],in_set);
      cmp = binary_seqs_cmp(outseq, outlist[i].plen, inseq, inlist[j].plen);
      if(cmp < 0) i++;
      else if(cmp == 0) {
        _update_store(&outlist[i], out_set, &inlist[j], in_set, store); j++;
      } else outlist[x++] = outlist[j++];
    }
  }

  for(; j < num_in; j++) {
    cset_in = path_set_colset(&inlist[j], in_set);
    if(!packedpath_is_colset_zero(cset_in, cbytes)) inlist[x++] = inlist[j];
  }

  in_set->members.len = x;
}

// Sum of bytes required to store set
size_t path_set_get_bytes_sum(const PathSet *set)
{
  size_t i, sum = 0;
  for(i = 0; i < set->members.len; i++) {
    sum += packedpath_mem2(set->cbytes, (set->members.data[i].plen+3)/4);
  }
  return sum;
}

// Load 'count' values into the set. Paths without the given colour are given
// a count of zero
void path_set_load_counts(PathSet *set, const PathStore *ps, size_t colour)
{
  ctx_assert(ps->colset_bytes == set->cbytes);
  ctx_assert2(ps->extra_bytes == 1, "PathStore doesn't hold counts");

  PathEntry *entry = set->members.data;
  const PathEntry *end = entry + set->members.len;
  const uint8_t *path;

  for(; entry < end; entry++)
  {
    if(bitset_get(path_set_colset(entry,set),colour)) {
      path = &ps->store[entry->pindex];
      entry->count = *(path + packedpath_mem(path,set->cbytes));
    }
    else entry->count = 0;
  }
}

// Get path with index before `last` which is in the colour `colour`
static inline int _get_prev_in_col(PathSet *set, int last, size_t colour)
{
  int i;
  for(i = last-1; i >= 0; i--) {
    if(bitset_get(path_set_colset(&set->members.data[i],set), colour))
      return i;
  }

  return -1;
}

// Remove colours from paths with 'count' < threshold
// Updates `pstore` if != NULL
void path_set_threshold(PathSet *set, uint8_t threshold,
                        size_t colour, uint8_t *pstore)
{
  if(threshold == 0) return;
  if(set->members.len == 0) return;
  if(!set->sorted) path_set_sort(set);

  int num_members = set->members.len;
  PathEntry *entries = set->members.data, *prev;
  const PathEntry *curr;

  // Work backwards sum up counts
  int j = _get_prev_in_col(set, num_members, colour);
  int i = _get_prev_in_col(set, j, colour);

  while(i >= 0)
  {
    prev = &entries[i];
    curr = &entries[j];

    if(prev->orient == curr->orient && prev->plen <= curr->plen &&
       binary_seqs_cmp(path_set_seq(prev, set), prev->plen,
                       path_set_seq(curr, set), prev->plen) == 0)
    {
      // path i is subset of j
      prev->count = MIN2(255, (size_t)prev->count + curr->count);
    }

    j = i;
    i = _get_prev_in_col(set, j, colour);
  }

  // Remove colours from low scoring paths in both this PathSet and PathStore
  for(i = 0; i < num_members; i++) {
    if(entries[i].count > 0 && entries[i].count < threshold) {
      bitset_del(path_set_colset(&entries[i], set), colour);
      if(pstore != NULL)
        bitset_del(packedpath_get_colset(pstore + entries[i].pindex), colour);
    }
  }
}

// Run checks on a PathSet
// Prints error message and returns false on error
// Otherwise returns true
bool path_set_check(const PathSet *set)
{
  size_t i;
  bool no_errors = true;
  for(i = 0; i < set->members.len; i++) {
    const PathEntry *entry = &set->members.data[i];
    if(entry->seq + (entry->plen+3)/4 > set->seqs.len) {
      warn("[PathSet] path %zu looks corrupted", i);
      no_errors = false;
    }
  }
  return no_errors;
}
