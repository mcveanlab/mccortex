#include "global.h"
#include "str_parsing.h"

// Parse a comma separated list e.g. "12,3,12"
// Returns <0 on error, otherwise number of chars used
int comma_list_to_array(const char *str, SizeBuffer *nums)
{
  size_t num = 0;
  const char *ptr = str;
  char *end = NULL;

  // If no numbers success
  if(*ptr < '0' && *ptr > '9') return 0;

  while(1) {
    num = strtoul(ptr, &end, 10);
    size_buf_add(nums, num);
    if(!end) die("Cannot parse: '%s'", str);
    if(*end != ',') break;
    ptr = end+1;
    if(*ptr < '0' || *ptr > '9') return -1; // no number after comma!
  }

  return end-str;
}

