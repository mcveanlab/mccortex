#include "global.h"
#include "util.h"

#include <math.h>

// comparison returns:
//   negative iff a < b
//          0 iff a == b
//   positive iff a > b

#define cmpfunc(fname,type_t)                                                  \
int fname(const void *a, const void *b) {                                      \
  const type_t *a2 = (const type_t *)a;                                        \
  const type_t *b2 = (const type_t *)b;                                        \
  return (*a2 < *b2 ? -1 : (*a2 > *b2));                                       \
}

cmpfunc(cmp_int, int);
cmpfunc(cmp_long, long);
cmpfunc(cmp_float, float);
cmpfunc(cmp_double, double);
cmpfunc(cmp_uint32, uint32_t);
cmpfunc(cmp_uint64, uint64_t);
cmpfunc(cmp_size, size_t);
cmpfunc(cmp_ptr, void*);


char parse_entire_int(char *str, int *result)
{
  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || *strtol_last_char_ptr == '\0')
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
  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(tmp <= UINT_MAX && *strtol_last_char_ptr == '\0')
  {
    *result = (unsigned int)tmp;
    return 1;
  }
  return 0;
}

char parse_entire_ulong(char *str, unsigned long *result)
{
  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(*strtol_last_char_ptr == '\0')
  {
    *result = tmp;
    return 1;
  }
  return 0;
}

// Load a string of numbers into an array. Separator can be any non-numerical
// character. e.g. '1,2,3' '1 2 3'
// Returns number of unsigned integers parsed.
// Sets *more = 1 if string end not reached
uint32_t parse_uint_liststr(const char *str, uint32_t *arr, uint32_t arrlen,
                            boolean *more)
{
  const char *tmp = str;
  char *endptr = NULL;
  uint32_t i;
  for(i = 0; i < arrlen; i++)
  {
    while(*tmp != '\0' && (*tmp < '0' || *tmp > '9')) tmp++;
    if(*tmp == '\0') break;
    arr[i] = strtol(tmp, &endptr, 10);
    tmp = endptr;
  }
  *more = (*tmp != '\0');
  return i;
}

uint32_t len_uint_liststr(const char *str)
{
  const char *tmp = str;
  uint32_t i;
  for(i = 0; ; i++)
  {
    while(*tmp != '\0' && (*tmp < '0' || *tmp > '9')) tmp++;
    if(*tmp == '\0') break;
    while(*tmp >= '0' && *tmp <= '9') tmp++;
  }
  return i;
}

boolean parse_uint_liststr_strict(const char *str, char sep,
                                  uint32_t *arr, uint32_t arrlen)
{
  const char *tmp = str;
  char *endptr = NULL;
  uint32_t i;
  for(i = 0; i < arrlen; i++)
  {
    if(*tmp < '0' || *tmp > '9') return false;
    arr[i] = strtol(tmp, &endptr, 10);
    tmp = endptr;
    if(i+1 < arrlen) {
      if(*tmp != sep) return false;
      tmp++;
    }
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
  if(strcasecmp(endptr,"T") == 0 || strcasecmp(endptr,"TB") == 0)
  { *bytes = num<<40; return true; }
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

/* Formating Numbers */

unsigned int num_of_digits(unsigned long num)
{
  unsigned int digits;
  for(digits = 1; num >= 10; digits++) num /= 10;
  return digits;
}

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27
// returns pointer to result
char* ulong_to_str(unsigned long num, char* result)
{
  int digits = num_of_digits(num);
  int i, num_commas = (digits-1) / 3;
  char *p = result + digits + num_commas;
  *(p--) = '\0';

  for(i = 0; i < digits; i++, num /= 10) {
    if(i > 0 && i % 3 == 0) *(p--) = ',';
    *(p--) = '0' + (num % 10);
  }

  return result;
}

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('-9,223,372,036,854,775,808')+1 = 27
char* long_to_str(long num, char* result)
{
  if(num < 0) {
    result[0] = '-';
    ulong_to_str(-num, result+1);
  }
  else
    ulong_to_str(num, result);

  return result;
}

// result must be long enough for result + 1 ('\0').
// Max length required is: 26+1+decimals+1 = 28+decimals bytes
// strlen('-9,223,372,036,854,775,808') = 27
// strlen('.') = 1
// +1 for \0
char* double_to_str(double num, int decimals, char* str)
{
  unsigned long whole_units = (unsigned long)num;
  num -= whole_units;

  ulong_to_str(whole_units, str);

  if(decimals > 0)
  {
    // Horrible hack to save character being overwritten with a leading zero
    // e.g. 12.121 written as '12' then '0.121', giving '10.121', put back '2'
    // '12.121'
    size_t offset = strlen(str);
    char c = str[offset-1];
    sprintf(str+offset-1, "%.*lf", decimals, num);
    str[offset-1] = c;
  }

  return str;
}

// str must be 26 + 3 + 1 + num decimals + 1 = 31+decimals bytes
// breakdown:
// strlen('18,446,744,073,709,551,615') = 26
// strlen(' GB') = 3
// strlen('.') = 1
// +1 for '\0'
char* bytes_to_str(unsigned long num, int decimals, char* str)
{
  const unsigned int num_unit_sizes = 7;
  const char *units[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB"};

  unsigned long unit;
  unsigned long num_cpy = num;

  for(unit = 0; num_cpy >= 1024 && unit < num_unit_sizes; unit++)
    num_cpy /= 1024;

  unsigned long bytes_in_unit = 0x1UL << (10 * unit);
  double num_of_units = (double)num / bytes_in_unit;

  double_to_str(num_of_units, decimals, str);
  char *ptr = str+strlen(str)-1;
  if(decimals > 0) {
    // Trim excess zeros
    while(ptr > str && *ptr == '0') ptr--;
    if(*ptr == '.') ptr--;
  }
  strcpy(ptr+1, units[unit]);

  return str;
}

// Number to string using G to mean 10^9, M to mean 10^6 etc
char* num_to_str(unsigned long num, int decimals, char* str)
{
  bytes_to_str(num, decimals, str);
  // Trim 'B' From the end
  size_t len = strlen(str);
  str[len-1] = '\0';

  return str;
}

/*
// This is needed if POSIX string functions not available
char* strdup(const char *str)
{
  size_t n = strlen(str);
  char *dup = malloc2(n+1);
  if(dup) memcpy(dup, str, (n+1)*sizeof(char));
  return dup;
}
*/

//
// Maths
//

// log(n!) = sum from i=1 to n, of (log(i))
float log_factorial(unsigned int number)
{
  unsigned int i;
  float ret = 0;
  for(i = 1; i <= number; i++) ret += log(i);
  return ret;
}

float log_factorial_ll(unsigned long long number)
{
  unsigned long long i;
  float ret = 0;
  for(i = 1; i <= number; i++) ret += log(i);
  return ret;
}

unsigned long calculate_mean_ulong(unsigned long *array, unsigned long len)
{
  unsigned long i, sum = 0, num = 0;

  for(i = 0; i < len; i++)
  {
    sum += i * array[i];
    num += array[i];
  }

  return num == 0 ? 0 : (sum / num);
}


//
// Time
//

// output of form: "10 days 23 hours 59 mins 59 secs"
// extreme: "-2147483647 days 23 hours 59 mins 59 secs",
// so str should be at least 42 bytes long
// returns number of bytes written
size_t seconds_to_str(unsigned long seconds, char *str)
{
  char *ptr = str;
  int days, hours, mins;
  days = seconds / (60 * 60 * 24);
  seconds -= days * (60 * 60 * 24);
  hours = seconds / (60 * 60);
  seconds -= hours * (60 * 60);
  mins = seconds / 60;
  seconds -= mins * 60;
  if(days > 0)
    ptr += sprintf(ptr, "%i day%s ", days, days == 1 ? "" : "s");
  if(days+hours > 0)
    ptr += sprintf(ptr, "%i hour%s ", hours, hours == 1 ? "" : "s");
  if(days+hours+mins > 0)
    ptr += sprintf(ptr, "%i min%s ", mins, mins == 1 ? "" : "s");
  ptr += sprintf(ptr, "%i secs", (int)seconds);
  return ptr - str;
}

//
// Genetics
//
const char complement_base[128] = {0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,'T',0,'G',  0,0,0,'C',0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,'A',0,0,  0,0,0,0,0,0,0,0,0,
                                   0,'t',0,'g',  0,0,0,'c',0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,'a',0,0,  0,0,0,0,0,0,0,0,0};

#ifndef NDEBUG
// These have asserts in them
char char_nucleotide_complement(char c)
{
  char cmplmnt = complement_base[(int)c];
  assert(cmplmnt != 0);
  return cmplmnt;
}
#endif

// length is the length in number of bases
// the char* should have one MORE base than that allocated, to hold '\0'
char *reverse_complement_str(char *str, size_t length)
{
  assert(strlen(str) >= length);

  if(length == 0) return str;
  if(length == 1) { str[0] = char_nucleotide_complement(str[0]); return str; }

  size_t i, j;
  for(i = 0, j = length - 1; i <= j; i++, j--)
  {
    char tmp = str[i];
    str[i] = char_nucleotide_complement(str[j]);
    str[j] = char_nucleotide_complement(tmp);
  }

  return str;
}
