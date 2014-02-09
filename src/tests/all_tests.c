#include "global.h"
#include "all_tests.h"

// Common functions here
FILE *ctx_tst_out = NULL;

size_t tests_num_run = 0, tests_num_failed = 0;

//
// Useful functions
//
void fill_rand(uint8_t *arr, size_t n)
{
  uint32_t r; size_t i;
  for(i = 0; i+4 <= n; i+=4) { r = rand(); memcpy(&arr[i], &r, 4); }
  r = rand();
  for(; i < n; i++) { arr[i] = (uint8_t)(r&0xff); r >>= 8; }
}

void rand_nucs(Nucleotide *nucs, size_t len)
{
  if(!len) return;
  size_t i, r = 0;
  for(i = 0; i < len; i++) {
    if((i & 15) == 0) r = (size_t)rand(); // 2 bits per cycle, 32 bits in rand()
    nucs[i] = r&3;
    r >>= 2;
  }
}

void bitarr_tostr(const uint8_t *arr, size_t len, char *str)
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
