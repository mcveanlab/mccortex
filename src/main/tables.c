#include <stdlib.h>
#include <stdio.h>

// Swap lowest four bits. A nibble is 4 bits (i.e. half a byte)
#define rev_nibble(x) ((((x)&1)<<3)|(((x)&2)<<1)|(((x)&4)>>1)|(((x)&8)>>3))

int main()
{
  size_t i;
  printf("static const uint8_t rev_nibble_arr[16]\n");
  printf("  = {");
  for(i = 0; i < 16; i++)
    printf("%s%zu", i ? ", " : "", rev_nibble(i));
  printf("}\n");
  return EXIT_SUCCESS;
}

