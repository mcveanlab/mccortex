#include "global.h"
#include "packed_path.h"
#include "all_tests.c"

static void fill_rand(uint8_t *arr, size_t n)
{
  int r; size_t i;
  for(i = 0; i+4 < n; i+=4) { r = rand(); memcpy(&arr[i], &r, 4); }
  r = rand();
  for(i = 4*(n/4); i<n; i++) { arr[i] = (uint8_t)(r&0xff); r >>= 8; }
}

static void bitarr_tostr(const uint8_t *arr, size_t len, char *str)
{
  size_t i, j;
  for(i = len-1; i != SIZE_MAX; i--) {
    for(j = 7; j != SIZE_MAX; j--) {
      *str = (arr[i]&(1U<<j)) ? '1' : '0';
      str++;
    }
  }
  *str = '\0';
}

void test_packed_path()
{
  status("[packedpath] Testing shift copy");

  uint8_t d0[10] = {0,0,0,0,0,0,0,0,0,0};
  uint8_t out[100];
  size_t i, j, shift, len;
  size_t t;

  // Shifting an array of zeros results in zeros
  // surrounding array should remain all ones
  memset(out, 0xff, 100);
  for(shift = 0; shift < 4; shift++) {
    packed_cpy_fast(out+1, d0, shift, 15); // first 4 bytes
    assert(out[0]==0xff);
    for(i = 1; i < 5; i++) assert(out[i]==0);
    for(i = 5; i < 100; i++) assert(out[i]==0xff);
  }

  #define NTESTS 100
  #define TLEN 100

  // Random testing
  uint8_t in[TLEN], slow[TLEN], med[TLEN], fast[TLEN];
  char tmp1[TLEN*8+1], tmp2[TLEN*8+1];

  for(t = 0; t < NTESTS; t++) {
    memset(slow, 0xff, TLEN);
    memset(med,  0xff, TLEN);
    memset(fast, 0xff, TLEN);
    fill_rand(in, TLEN);

    len = rand() % (TLEN/2+1);
    shift = rand() % 4;
    // printf("len: %zu shift: %zu\n", len, shift);

    packed_cpy_slow(slow, in, shift, len);
    packed_cpy_med(med, in, shift, len);
    packed_cpy_fast(fast, in, shift, len);

    if(len > shift) {
      // Check with string method to be extra safe
      bitarr_tostr(in,   TLEN, tmp1);
      bitarr_tostr(slow, TLEN, tmp2);
      for(i = 8*TLEN-1; i != 8*TLEN-(len-shift)*2; i--)
        assert(tmp1[i-shift*2] == tmp2[i]);
    }

    // Check all results match
    for(i = 0; i < TLEN && slow[i] == med[i]; i++);
    for(j = 0; j < TLEN && med[j] == fast[j]; j++);

    // Print output if arrays don't match
    if(i < TLEN || j < TLEN) {
      printf("len: %zu shift: %zu\n", len, shift);
      bitarr_tostr(in,   TLEN, tmp1); printf("in:  %s\n", tmp1);
      bitarr_tostr(slow, TLEN, tmp1); printf("slw: %s\n", tmp1);
      bitarr_tostr(med,  TLEN, tmp1); printf("med: %s\n", tmp1);
      bitarr_tostr(fast, TLEN, tmp1); printf("fst: %s\n", tmp1);
      printf("\n");
    }

    assert(i == TLEN);
    assert(j == TLEN);
  }
}
