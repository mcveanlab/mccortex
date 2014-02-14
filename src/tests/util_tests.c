#include "global.h"
#include "all_tests.h"
#include "util.h"

#include <math.h> // NAN, INFINITY

static void test_util_num_to_str()
{
  char str[100];
  // Check NaN and Inf are correctly written
  TASSERT2(strcmp(num_to_str(NAN, 2, str),"NaN") == 0, "Got: %s", str);
  TASSERT2(strcmp(num_to_str(INFINITY, 2, str),"Inf") == 0, "Got: %s", str);
}

static void test_util_bytes_to_str()
{
  test_status("[util] testing bytes_to_str()");

  char str[100];

  // Excess decimal points are trimmed off
  // 14.0MB -> 14MB
  TASSERT2(strcmp(bytes_to_str(14688256,1,str),"14MB") == 0, "Got: %s", str);
  // 1.9GB -> 1.9GB
  TASSERT2(strcmp(bytes_to_str(2040110000,1,str),"1.9GB") == 0, "Got: %s", str);
  // 1.99GB -> 2GB
  TASSERT2(strcmp(bytes_to_str(2140110000,1,str),"2GB") == 0, "Got: %s", str);
  // 1500KB -> 1.4MB
  TASSERT2(strcmp(bytes_to_str(1500000,1,str),"1.4MB") == 0, "Got: %s", str);
  // 0.5GB -> 512MB
  TASSERT2(strcmp(bytes_to_str(536900000,1,str),"512MB") == 0, "Got: %s", str);
  // 1 -> 1B
  TASSERT2(strcmp(bytes_to_str(1,1,str),"1B") == 0, "Got: %s", str);
  // 1023 -> 1023B
  TASSERT2(strcmp(bytes_to_str(1023,1,str),"1,023B") == 0, "Got: %s", str);
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
  test_util_num_to_str();
  test_util_bytes_to_str();
  test_util_get_GCD();
}
