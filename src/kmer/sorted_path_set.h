#ifndef SORTED_PATH_SET_H_
#define SORTED_PATH_SET_H_

#include "packed_path.h"
#include "path_store.h"
#include "hash_table.h" // defines dBNode

typedef struct
{
  uint8_t *seq; // colset is stored at (seq-cbytes)
  PathIndex pindex;
  PathLen orient:1, plen:PATH_LEN_BITS;
} PathEntry;

#include "objbuf_macro.h"
create_objbuf(pentrybuf, PathEntryBuffer, PathEntry);

#ifndef BYTE_BUFFER_DEFINED
  create_objbuf(bytebuf, ByteBuffer, uint8_t);
  #define BYTE_BUFFER_DEFINED
#endif

typedef struct
{
  hkey_t hkey;
  size_t cbytes; // number of bytes in colourset
  ByteBuffer seqs;
  PathEntryBuffer members;
} SortedPathSet;

#define sorted_path_colset(entry,set) ((entry)->seq - (set)->cbytes)

void sorted_path_set_alloc(SortedPathSet *set);
void sorted_path_set_dealloc(SortedPathSet *set);

void sorted_path_set_init(SortedPathSet *set, const PathStore *ps, hkey_t hkey);

void sorted_path_set_init2(SortedPathSet *set, hkey_t hkey, size_t cbytes,
                           const uint8_t *store, PathIndex pindex,
                           const FileFilter *fltr);

// Remove redundant entries
void sorted_path_set_slim(SortedPathSet *set);

// Remove entries from set that are in the filter set
// if `rmsubstr` then also remove substring matches
// Note: does not remove duplicates from `set`,
//       call sorted_path_set_slim() to do that
// if passed store, we use it to update colourset of paths in out_set
void sorted_path_set_merge(SortedPathSet *out_set, SortedPathSet *in_set,
                           bool rmsubstr, uint8_t *store);

#endif /* SORTED_PATH_SET_H_ */
