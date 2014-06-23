#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <limits.h>

/*
  Generate static lookup tables used in ninja-cortex
*/

// Swap lowest four bits. A nibble is 4 bits (i.e. half a byte)
#define rev_nibble(x) ((((x)&1)<<3)|(((x)&2)<<1)|(((x)&4)>>1)|(((x)&8)>>3))

static void rev_nibble_table()
{
  size_t i;
  printf("const uint8_t rev_nibble_table[16]\n");
  printf("  = {");
  for(i = 0; i < 16; i++)
    printf("%s%zu", i ? ", " : "", rev_nibble(i));
  printf("};\n\n");
}

static void rev_cmp_table()
{
  printf("const uint8_t revcmp_table[256] = \n{\n");
  uint8_t a, b, c, d, x, col = 0, i = 0;

  for(a = 0; a < 4; a++) {
    for(b = 0; b < 4; b++) {
      for(c = 0; c < 4; c++) {
        for(d = 0; d < 4; d++) {
          x = ~((d<<6)|(c<<4)|(b<<2)|a);
          printf("%s", col ? " " : "  ");
          printf("0x%02zx", (size_t)x);
          if(i < 255) printf(",");
          if(col == 7) { printf("\n"); col = 0; }
          else col++;
          i++;
        }
      }
    }
  }
  printf("};\n\n");
}

// static const uint8_t nibble_popcount_table[16] = {0,1,1,2,1,2,2,3,1,12,2,3,2,3,3,4};
static void count_bits_table()
{
  int i;
  printf("static const uint8_t nibble_popcount_table[16] = {0");
  for(i = 1; i < 16; i++)
    printf(",%i", __builtin_popcount(i));
  printf("};\n\n");
}

int main()
{
  rev_nibble_table();
  rev_cmp_table();
  count_bits_table();
  return EXIT_SUCCESS;
}

