#include "global.h"
#include "all_tests.h"
#include "util.h"

static void test_util_bytes_to_str()
{
  test_status("[util] testing bytes_to_str()");

  char str[100];
  // Excess decimal points are trimmed off 14.0MB -> 14MB
  assert(strcmp(bytes_to_str(14688256,1,str),"14MB") == 0);
}

static void test_util_get_GCD()
{
  test_status("[util] testing get_GCD()");
  assert(calc_GCD(0,0) == 0);
  assert(calc_GCD(10,0) == 10);
  assert(calc_GCD(0,10) == 10);
  assert(calc_GCD(2,2) == 2);
  assert(calc_GCD(1,1) == 1);
  assert(calc_GCD(1,2) == 1);
  assert(calc_GCD(1,100) == 1);
  assert(calc_GCD(2,4) == 2);
  assert(calc_GCD(7,5) == 1);
  assert(calc_GCD(18,6) == 6);
  assert(calc_GCD(3,6) == 3);
  assert(calc_GCD(100,120) == 20);
  assert(calc_GCD(100,125) == 25);
}

void test_util()
{
  test_util_bytes_to_str();
  test_util_get_GCD();
}
