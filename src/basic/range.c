#include "global.h"
#include "range.h"

// If range is empty assume 0..(N-1) inclusive

// single int or range e.g. '6' [ans: 1] or '2-12' [ans: 10], '*' means all
// ranges are allowed to go backwards e.g. 10-8 means 10,9,8
// returns number of bytes consumed or -1 if error
static int range_parse(const char *range_str, size_t *start, size_t *end,
                       size_t range_max)
{
  char *endptr;
  const char *str = range_str;
  unsigned long from, to;

  if(*str == '*') {
    *start = 0;
    *end = range_max;
    return 1;
  }

  from = strtoul(str, &endptr, 10);
  to = from;
  if(endptr == str) return -1;
  if(*endptr == '-') {
    str = endptr+1;
    to = strtoul(str, &endptr, 10);
    if(endptr == str) return -1;
  }
  if(from > range_max || to > range_max) return -1;
  *start = from;
  *end = to;

  return endptr - range_str;
}

size_t range_get_num(const char *str, size_t range_max)
{
  const char *ptr = str;
  size_t start, end, num_cols = 0;
  int bytes;

  while(*ptr != '\0')
  {
    if((bytes = range_parse(ptr, &start, &end, range_max)) == -1)
      die("Invalid range specifier: %s [max: %zu]", str, range_max);
    ptr += bytes;
    num_cols += ABSDIFF(start,end)+1;
    if(*ptr == ',') ptr++;
  }

  return num_cols == 0 ? range_max+1 : num_cols;
}

void range_parse_array(const char *str, size_t *arr, size_t range_max)
{
  const char *ptr = str;
  size_t num_cols = 0, j, start, end;
  int bytes;

  while(*ptr != '\0')
  {
    if((bytes = range_parse(ptr, &start, &end, range_max)) == -1)
      die("Invalid range specifier: %s [max: %zu]", str, range_max);
    ptr += bytes;
    if(*ptr == ',') ptr++;
    if(start <= end)
      for(j = start; j <= end; j++) arr[num_cols++] = j;
    else
      for(j = start; j <= start; j--) arr[num_cols++] = j;
  }

  if(ptr > str && *(ptr-1) == ',')
      die("Invalid range specifier: %s [max: %zu]", str, range_max);

  if(num_cols == 0) {
    for(num_cols = 0; num_cols <= range_max; num_cols++)
      arr[num_cols] = num_cols;
  }
}
