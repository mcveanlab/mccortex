#include "global.h"
#include "all_tests.h"
#include "util.h"

static void test_util_bytes_to_str()
{
  test_status("[util] testing bytes_to_str()");

  char str[100];

  // Excess decimal points are trimmed off 14.0MB -> 14MB
  bytes_to_str(14688256,1,str);
  TASSERT2(strcmp(str,"14MB") == 0, "Got: %s", str);
}

static void test_util_get_GCD()
{
  test_status("[util] testing get_GCD()");
  TASSERT(calc_GCD(0,0) == 0);
  TASSERT(calc_GCD(10,0) == 10);
  TASSERT(calc_GCD(0,10) == 10);
  TASSERT(calc_GCD(2,2) == 2);
  TASSERT(calc_GCD(1,1) == 1);
  TASSERT(calc_GCD(1,2) == 1);
  TASSERT(calc_GCD(1,100) == 1);
  TASSERT(calc_GCD(2,4) == 2);
  TASSERT(calc_GCD(7,5) == 1);
  TASSERT(calc_GCD(18,6) == 6);
  TASSERT(calc_GCD(3,6) == 3);
  TASSERT(calc_GCD(100,120) == 20);
  TASSERT(calc_GCD(100,125) == 25);
}

void test_util()
{
  test_util_bytes_to_str();
  test_util_get_GCD();
}
