#include "global.h"
#include "all_tests.h"
#include "packed_path.h"
#include "binary_seq.h"

#define NTESTS 100

void test_len_orient()
{
  test_status("[packedpath] Testing combine_lenorient");

  PathLen len, len2, merged;
  Orientation orient, orient2;
  int r; size_t t;

  for(t = 0; t < NTESTS; t++) {
    r = rand();
    len = r & PP_LENMASK;
    orient = r >> 31;
    merged = packedpath_combine_lenorient(len,orient);
    len2 = packedpath_len(merged);
    orient2 = packedpath_or(merged);
    TASSERT(len == len2);
    TASSERT(orient == orient2);
  }
}

void test_packed_path()
{
  test_len_orient();
}
