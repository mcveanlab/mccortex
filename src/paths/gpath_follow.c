#include "global.h"
#include "gpath_follow.h"
// #include "binary_seq.h"

/*
// Check if the GPathFollow cache needs updated, based off path->pos value
// if it does, update it
void gpath_follow_cache_update(GPathFollow *path, size_t pos)
{
  size_t fetch_offset, fetch_bytes, total_bytes;

  // 4 bases per byte
  #define BASES_IN_CACHE(p) (sizeof((p)->cache) * 4)
  uint16_t new_cache_start = (pos / BASES_IN_CACHE(path)) * BASES_IN_CACHE(path);

  if(new_cache_start != path->first_cached)
  {
    path->first_cached = new_cache_start;
    fetch_offset = path->first_cached/4;
    total_bytes = binary_seq_mem(path->len);
    fetch_bytes = MIN2(total_bytes-fetch_offset, sizeof(path->cache));
    memcpy(path->cache, path->gpath->seq + fetch_offset, fetch_bytes);
    memset(path->cache+fetch_bytes, 0, sizeof(path->cache)-fetch_bytes);
    // Need to zero rest of cache since it is used in hashing
    //  -> must be deterministic
  }
}

// Get a base from the GPathFollow cache
Nucleotide gpath_follow_get_base(GPathFollow *path, size_t pos)
{
  ctx_assert2(pos < path->len, "pos: %i path->len: %i", (int)pos, (int)path->len);
  gpath_follow_cache_update(path, pos);
  return binary_seq_get(path->cache, pos - path->first_cached);
}
*/

// For GraphWalker to work we assume all edges are merged into one colour
// (i.e. graph->num_edge_cols == 1)
// If only one colour loaded we assume all edges belong to this colour

GPathFollow gpath_follow_create(const GPath *gpath)
{
  GPathFollow fpath = {.gpath = gpath,
                        .pos = 0,
                       .len = gpath->num_juncs,
                       .age = 0};

  // .first_cached = 1 is invalid (not multiple of sizeof(cache)*4), so forces
  // fetch on first request
  // fpath.first_cached = 1;
  // memset(fpath.cache, 0, sizeof(fpath.cache));

  return fpath;
}
