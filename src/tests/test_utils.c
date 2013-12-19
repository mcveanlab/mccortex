#include "global.h"
#include "all_tests.h"

void rnd_seq(char *seq, size_t len)
{
  const char bases[4] = "ACGT";
  size_t i, r;

  for(i = 0; i < len; i++) {
    if((i & 15) == 0) r = rand();
    seq[i] = bases[r & 3]; r >>= 2;
  }

  seq[len] = '\0';
}
