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

static void rand_nucs(Nucleotide *nucs, size_t len)
{
  if(!len) return;
  size_t i, r;
  for(i = 0; i < len; i++) {
    if((i & 15) == 0) r = (size_t)rand(); // 2 bits per cycle, 32 bits in rand()
    nucs[i] = r&3;
    r >>= 2;
  }
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

#define NTESTS 100
#define TLEN 200

static void test_pack_cpy()
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

static void print_nucs(Nucleotide *nucs, size_t len) {
  size_t i;
  for(i = 0; i < len; i++) printf(" %u", (uint32_t)nucs[i]);
}

static void test_pack_unpack()
{
  status("[packedpath] Testing pack_bases() / unpack_bases()");

  uint8_t packed[TLEN];
  Nucleotide bases0[TLEN], bases1[TLEN];
  char str[8*TLEN+1];
  size_t i, t, len;

  // Run NTESTS
  // randomize bases0, pack into packed, unpack into bases1
  // compare bases0 vs bases1
  for(t = 0; t < NTESTS; t++) {
    len = rand() % (TLEN/2+1);
    rand_nucs(bases0, len);
    memset(packed, 0, TLEN);
    pack_bases(packed, bases0, len);
    unpack_bases(packed, bases1, len);

    for(i = 0; i < len && bases0[i] == bases1[i]; i++);

    // print output if input != output
    if(i != len) {
      bitarr_tostr(packed, (len*2+7)/8, str);
      printf("bases0: "); print_nucs(bases0, len); printf("\n");
      printf("bases1: "); print_nucs(bases1, len); printf("\n");
      printf("packed: %s\n", str);
    }

    assert(i == len);
  }
}

void test_packed_path()
{
  test_pack_cpy();
  test_pack_unpack();
}
