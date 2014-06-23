#include "global.h"
#include "gpath_follow.h"
#include "binary_seq.h"

// Check if the GPathFollow cache needs updated, based of path->pos value
// if it does, update it
void gpath_follow_cache_update(GPathFollow *path, size_t pos)
{
  size_t fetch_offset, fetch_bytes, total_bytes;

  // 4 bases per byte
  uint16_t new_cache_start = sizeof(path->cache) * 4 *
                             (pos / (sizeof(path->cache) * 4));

  if(new_cache_start != path->first_cached)
  {
    path->first_cached = new_cache_start;
    fetch_offset = path->first_cached/4;
    total_bytes = (path->len+3)/4;
    fetch_bytes = MIN2(total_bytes-fetch_offset, sizeof(path->cache));
    memcpy(path->cache, path->gpath->seq + fetch_offset, fetch_bytes);
    memset(path->cache+fetch_bytes, 0, sizeof(path->cache)-fetch_bytes);
  }
}

// Get a base from the GPathFollow cache
Nucleotide gpath_follow_get_base(GPathFollow *path, size_t pos)
{
  gpath_follow_cache_update(path, pos);
  return binary_seq_get(path->cache, pos - path->first_cached);
}

// For GraphWalker to work we assume all edges are merged into one colour
// (i.e. graph->num_edge_cols == 1)
// If only one colour loaded we assume all edges belong to this colour

GPathFollow gpath_follow_create(const GPath *gpath)
{
  // .first_cached = 1 is invalid (not multiple of sizeof(cache)*4), so forces
  // fetch on first request
  GPathFollow fpath = {.gpath = gpath, .pos = 0, .len = gpath->num_juncs,
                       .first_cached = 1};
  memset(fpath.cache, 0, sizeof(fpath.cache));
  return fpath;
}
