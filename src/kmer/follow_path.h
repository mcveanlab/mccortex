#ifndef FOLLOW_PATH_H_
#define FOLLOW_PATH_H_

#include "binary_seq.h"
#include "packed_path.h"

// This struct is packed so we can hash it quickly
struct FollowPathStruct
{
  const uint8_t *seq;
  PathLen pos, len;
  // A small buffer of upcoming 24 bases
  PathLen first_cached; // first base in buffer (multiple of 4: 0,4,8,...)
  uint8_t cache[6]; // first..first+24-1 (24 bases)
} __attribute__((packed));

typedef struct FollowPathStruct FollowPath;

#include "objbuf_macro.h"
create_objbuf(path_buf,PathBuffer,FollowPath);

// Check if the FollowPath cache needs updated, based of path->pos value
// if it does, update it
static inline void fpath_cache_update(FollowPath *path, size_t pos)
{
  size_t fetch_offset, fetch_bytes, total_bytes;

  // 4 bases per byte
  PathLen new_cache_start = sizeof(path->cache) * 4 *
                            (pos / (sizeof(path->cache) * 4));

  if(new_cache_start != path->first_cached)
  {
    path->first_cached = new_cache_start;
    fetch_offset = path->first_cached/4;
    total_bytes = (path->len+3)/4;
    fetch_bytes = MIN2(total_bytes-fetch_offset, sizeof(path->cache));
    memcpy(path->cache, path->seq + fetch_offset, fetch_bytes);
    memset(path->cache+fetch_bytes, 0, sizeof(path->cache)-fetch_bytes);
  }
}

// Get a base from the FollowPath cache
static inline Nucleotide fpath_get_base(FollowPath *path, size_t pos)
{
  fpath_cache_update(path, pos);
  return binary_seq_get(path->cache, pos - path->first_cached);
}

// For GraphWalker to work we assume all edges are merged into one colour
// (i.e. graph->num_edge_cols == 1)
// If only one colour loaded we assume all edges belong to this colour

static inline FollowPath fpath_create(const uint8_t *seq, PathLen plen)
{
  // .first_cached = 1 is invalid (not multiple of sizeof(cache)*4), so forces
  // fetch on first request
  FollowPath fpath = {.seq = seq, .pos = 0, .len = plen, .first_cached = 1};
  memset(fpath.cache, 0, sizeof(fpath.cache));
  return fpath;
}

#endif /* FOLLOW_PATH_H_ */
