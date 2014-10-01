#ifndef GPATH_FOLLOW_H_
#define GPATH_FOLLOW_H_

#include "dna.h"
#include "gpath.h"

// This struct is packed so we can hash it quickly
struct GPathFollowStruct
{
  const GPath *gpath;
  uint16_t pos, len;
  uint32_t age;
  // A small buffer of upcoming 24 bases
  uint16_t first_cached; // first base in buffer (multiple of 4: 0,4,8,...)
  uint8_t cache[6]; // first..first+23 (24 bases)
} __attribute__((packed));

typedef struct GPathFollowStruct GPathFollow;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(gpath_follow_buf,GPathFollowBuffer,GPathFollow);

void gpath_follow_cache_update(GPathFollow *path, size_t pos);
Nucleotide gpath_follow_get_base(GPathFollow *path, size_t pos);
GPathFollow gpath_follow_create(const GPath *gpath);

#endif /* GPATH_FOLLOW_H_ */
