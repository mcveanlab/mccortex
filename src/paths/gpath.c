#include "global.h"
#include "gpath.h"
#include "util.h"
#include "binary_seq.h"

// Compare by orient, sequence
int gpath_cmp(const GPath *a, const GPath *b)
{
  int ret;
  if((ret = (int)a->orient - (int)b->orient) != 0) return ret;
  if((ret = binary_seqs_cmp(a->seq, a->num_juncs, b->seq, b->num_juncs)) != 0)
    return ret;
  return 0;
}

size_t gpath_colset_bits_set(const GPath *gpath, size_t ncols)
{
  size_t i, nbytes = (ncols+7)/8, nbits_set = 0;
  const uint8_t *colset = gpath_get_colset(gpath, ncols);
  for(i = 0; i < nbytes; i++) nbits_set += byte_popcount(colset[i]);
  return nbits_set;
}

// Copy colour set bits from src to dst
void gpath_colset_or_mt(GPath *dst_gp, const GPath *src_gp, size_t ncols)
{
  uint8_t *dst = gpath_get_colset(dst_gp, ncols);
  const uint8_t *src = gpath_get_colset(src_gp, ncols);
  size_t i, nbytes = (ncols+7)/8;
  for(i = 0; i < nbytes; i++)
    __sync_fetch_and_or((volatile uint8_t*)&dst[i], src[i]);
}

// Remove from `set0` bits that are set in `set1`
// Returns 0 if no colours remain in src path, 1 otherwise
uint8_t gpath_colset_rm_intersect(const GPath *dst_gp, GPath *src_gp,
                                  size_t ncols)
{
  const uint8_t *dst = gpath_get_colset(dst_gp, ncols);
  uint8_t *src = gpath_get_colset(src_gp, ncols);
  const uint8_t *end; uint8_t src_or = 0;
  for(end = dst + (ncols+7)/8; dst < end; dst++, src++) {
    *src &= ~*dst;
    src_or |= *src;
  }
  return src_or ? 1 : 0;
}
