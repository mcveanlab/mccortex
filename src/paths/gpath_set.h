#ifndef GPATH_SET_H_
#define GPATH_SET_H_

#include "gpath.h"
#include "common_buffers.h"
#include "binary_seq.h"

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(gpath_buf, GPathBuffer, GPath);

typedef uint64_t pkey_t;

// These passed around to be added
typedef struct
{
  uint8_t *seq;
  uint8_t *colset, *nseen; // if not null, colset: ncol bits; nseen: ncol bytes
  uint32_t klen;
  uint16_t num_juncs;
  Orientation orient;
} GPathNew;

typedef struct
{
  const size_t ncols;
  GPathBuffer entries;
  ByteBuffer seqs; // colset+seq for each path
  ByteBuffer nseen_buf; // counts for how many times we've seen path
  Uint32Buffer klen_buf; // kmer length of buffer
  bool can_resize;
} GPathSet;


static inline pkey_t gpset_get_pkey(const GPathSet *gpset, const GPath *gpath)
{
  ctx_assert2(gpath >= gpset->entries.data &&
              gpath <= gpset->entries.data+gpset->entries.len,
              "GPath is not in GPathSet");
  return gpath - gpset->entries.data;
}

// If resize true, cannot do multithreaded but can resize array
// If resize false, die if out of mem, but can multithread
void gpath_set_alloc2(GPathSet *gpset, size_t ncols,
                      size_t initpaths, size_t initmem,
                      bool resize, bool keep_path_counts);

// If resize true, cannot do multithreaded but can resize array
// If resize false, die if out of mem, but can multithread
void gpath_set_alloc(GPathSet *set, size_t ncols, size_t initmem,
                     bool resize, bool keep_path_counts);
void gpath_set_dealloc(GPathSet *set);
void gpath_set_reset(GPathSet *set);

void gpath_set_print_stats(const GPathSet *gpset);

// Always adds new path. If newpath could be a duplicate, use gpathhash
// Threadsafe only if resize is false. GPath* not safe to edit until it returns
// Copies newgpath.seq over and wipe new colset
GPath* gpath_set_add_mt(GPathSet *gpset, GPathNew newgpath);

// Returns true if we are storing number of sightings and kmer length
#define gpath_set_has_nseen(gpset) ((gpset)->nseen_buf.data != NULL)

// Get kmer length of a GPath
uint32_t gpath_set_get_klen(const GPathSet *gpset, const GPath *gpath);

uint8_t* gpath_set_get_nseen(const GPathSet *gpset, const GPath *gpath);

// Copy nseen counts to dst from src
void gpath_set_nseen_sum_mt(const GPath *dst, GPathSet *dstset,
                            const GPath *src, const GPathSet *srcset);
// Copy nseen counts to dst from src
void gpath_set_nseen_sum2_mt(GPath *dst, GPathSet *dstset,
                             const uint8_t *src_nseen);


GPathNew gpath_set_get(const GPathSet *gpset, const GPath *gpath);

#define gpaths_are_equal(a,b) \
  ((a).orient == (b).orient && \
   binary_seqs_cmp((a).seq, (a).num_juncs, (b).seq, (b).num_juncs) == 0)

#endif /* GPATH_SET_H_ */
