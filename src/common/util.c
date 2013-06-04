
#include "global.h"
#include "util.h"

#include <stddef.h>
#include <math.h>

// integer comparison: returns:
//   negative iff a < b
//          0 iff a == b
//   positive iff a > b
int int_cmp(const void *a, const void *b)
{
  // casting pointer types
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;

  return (*ia  - *ib);
}

// Sort an array of addresses in ascending order (1,2,10,20)
// NULL values will be put at the front
//   negative iff a < b
//          0 iff a == b
//   positive iff a > b
int pointer_address_cmp(const void *a, const void *b)
{
  // ptrdiff_t is defined in stddef.h and is the datatype for differences
  // between pointers (which are defined with size size_t)
  ptrdiff_t cmp = *(const ptrdiff_t *)a - *(const ptrdiff_t *)b;

  // Convert to an int (don't want integer to overflow)
  if(cmp == 0) return 0;
  else return (cmp > 0 ? 1 : -1);
}

char parse_entire_int(char *str, int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (int)tmp;
    return 1;
  }
}

char parse_entire_uint(char *str, unsigned int *result)
{
  size_t len = strlen(str);
  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(tmp <= UINT_MAX && strtol_last_char_ptr == str+len)
  {
    *result = (unsigned int)tmp;
    return 1;
  }
  return 0;
}

char parse_entire_ulong(char *str, unsigned long *result)
{
  size_t len = strlen(str);
  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(strtol_last_char_ptr == str+len)
  {
    *result = tmp;
    return 1;
  }
  return 0;
}

// Load a string of numbers into an array. Separator can be any non-numerical
// character. e.g. '1,2,3' '1 2 3'
// Returns 1 if the number of uint parsed == explen, 0 otherwise
boolean parse_uint_liststr(const char *str, uint32_t *arr, uint32_t explen)
{
  const char *tmp = str;
  char *endptr = NULL;
  uint32_t i;
  for(i = 0; i < explen; i++)
  {
    if(*tmp < '0' || *tmp > '9') tmp++;
    long j = strtol(tmp, &endptr, 10);
    if(tmp == endptr || j < 0) return 0;
    arr[i] = j;
    tmp = endptr;
  }
  return (*tmp == '\0');
}

size_t count_char(const char *str, char c)
{
  const char *tmp = str;
  size_t count = 0;
  while((tmp = strchr(tmp, c)) != NULL) { count++; tmp++; }
  return count;
}

boolean bases_to_integer(const char *arg, size_t *bases)
{
  char *endptr;
  double num = strtod(arg, &endptr);
  if(endptr == arg) return false;
  if(strcasecmp(endptr,"G") == 0 || strcasecmp(endptr,"GB") == 0 ||
     strcasecmp(endptr,"Gbp") == 0) { *bases = num * 1000000000; return true; }
  if(strcasecmp(endptr,"M") == 0 || strcasecmp(endptr,"MB") == 0 ||
     strcasecmp(endptr,"Mbp") == 0) { *bases = num * 1000000; return true; }
  if(strcasecmp(endptr,"K") == 0 || strcasecmp(endptr,"KB") == 0 ||
     strcasecmp(endptr,"Kbp") == 0) { *bases = num * 1000; return true; }
  if(strcasecmp(endptr,"b") == 0 || *endptr == '\0') { *bases = num; return true; }
  return false;
}

boolean mem_to_integer(const char *arg, size_t *bytes)
{
  char *endptr;
  unsigned long num = strtoul(arg, &endptr, 10);
  if(endptr == arg) return false;
  if(strcasecmp(endptr,"G") == 0 || strcasecmp(endptr,"GB") == 0)
  { *bytes = num<<30; return true; }
  if(strcasecmp(endptr,"M") == 0 || strcasecmp(endptr,"MB") == 0)
  { *bytes = num<<20; return true; }
  if(strcasecmp(endptr,"K") == 0 || strcasecmp(endptr,"KB") == 0)
  { *bytes = num<<10; return true; }
  if(*endptr != '\0') return false;
  *bytes = num;
  return true;
}

/*
// This is needed if POSIX string functions not available
char* strdup(const char *str)
{
  size_t n = strlen(str);
  char *dup = malloc(n+1);
  if(dup) memcpy(dup, str, (n+1)*sizeof(char));
  return dup;
}
*/

//
// Maths
//

// log(n!) = sum from i=1 to n, of (log(i))
float log_factorial(int number)
{
  assert(number >= 0);

  int i;
  float ret = 0;
  for(i = 1; i <= number; i++) ret += log(i);

  return ret;
}

float log_factorial_ll(long long number)
{
  assert(number >= 0);

  long long i;
  float ret = 0;
  for(i = 1; i <= number; i++) ret += log(i);

  return ret;
}

unsigned long calculate_mean_ulong(unsigned long *array, unsigned long len)
{
  unsigned long sum = 0;
  unsigned long num = 0;
  unsigned long i;

  for(i = 0; i < len; i++)
  {
    sum += i * array[i];
    num += array[i];
  }

  return num == 0 ? 0 : (sum / num);
}
