#include "global.h"
#include "all_tests.h"
#include "util.h"

void test_util()
{
  status("[util] testing bytes_to_str()");

  char str[100];
  // Excess decimal points are trimmed off 14.0MB -> 14MB
  assert(strcmp(bytes_to_str(14688256,1,str),"14MB") == 0);
}
