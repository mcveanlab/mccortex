
#include "global.h"
#include "delta_array.h"

// 100,6,20,7,52,8
//  100,6 => 100 long, arr[0] = 6
//  20,7  => arr[20] = 7
//  52,8  => arr[52] = 8

// 0,0 for list of length zero

void delta_arr_alloc(delta_array_t *list)
{
  list->changes_capacity = 256;
  list->indices = malloc(list->changes_capacity * sizeof(uint32_t));
  list->values = malloc(list->changes_capacity * sizeof(uint32_t));
  list->len = list->num_changes = 0;
  list->arr_capacity = 1024;
  list->arr = malloc(list->arr_capacity * sizeof(uint32_t));
}

void delta_arr_dealloc(delta_array_t *list)
{
  free(list->indices);
  free(list->values);
  free(list->arr);
}

void delta_arr_reverse(delta_array_t *list)
{
  uint32_t i, j, tmp, len = list->len;
  for(i = 0; i < list->num_changes; i++)
  {
    j = list->num_changes - i - 1;
    list->indices[i] = len - list->indices[i] - 1;
    list->indices[j] = len - list->indices[j] - 1;
    SWAP(list->indices[i], list->indices[j], tmp);
    SWAP(list->values[i], list->values[j], tmp);
  }
  
  uint32_t prev = list->indices[0];
  list->indices[0] = 0;

  for(i = 1; i < list->num_changes; i++)
  {
    tmp = list->indices[i];
    list->indices[i] = prev+1;
    prev = tmp;
  }
}

void delta_array_unpack(delta_array_t *list)
{
  if(list->len > list->arr_capacity)
  {
    list->arr_capacity = ROUNDUP2POW(list->len);
    list->arr = realloc(list->arr, list->arr_capacity * sizeof(uint32_t));
  }

  uint32_t i, idx = 0, value = list->values[idx];
  uint32_t nxtidx = idx+1 < list->num_changes ? list->indices[idx+1] : list->len;

  for(i = 0; i < list->len; i++)
  {
    if(i == nxtidx) {
      idx++;
      value = list->values[idx];
      nxtidx = idx+1 < list->num_changes ? list->indices[idx+1] : list->len;
    }
    list->arr[i] = value;
  }
}

void delta_arr_print(delta_array_t *list, FILE *fh)
{
  fprintf(fh, "%u,%u", list->len, list->values[0]);
  size_t i;
  for(i = 1; i < list->num_changes; i++)
    fprintf(fh, ",%u,%u", list->indices[i], list->values[i]);
}

void delta_arr_gzprint(delta_array_t *list, gzFile gz)
{
  gzprintf(gz, "%u,%u", list->len, list->values[0]);
  size_t i;
  for(i = 1; i < list->num_changes; i++)
    gzprintf(gz, ",%u,%u", list->indices[i], list->values[i]);
}

void delta_arr_from_str(const char *str, delta_array_t *list)
{
  char *endptr;

  assert(str[0] >= '0' && str[0] <= '9');
  list->len = strtol(str, &endptr, 10);
  assert(endptr[0] == ',');
  endptr++;

  assert(endptr[0] >= '0' && endptr[0] <= '9');
  list->indices[0] = 0;
  list->values[0] = strtol(endptr, &endptr, 10);
  list->num_changes = 1;

  while(*endptr == ',')
  {
    if(list->num_changes == list->changes_capacity)
    {
      list->changes_capacity *= 2;
      size_t mem_uint = list->changes_capacity * sizeof(uint32_t);
      list->indices = realloc(list->indices, mem_uint);
      list->values = realloc(list->values, mem_uint);
    }

    endptr++;
    assert(endptr[0] >= '0' && endptr[0] <= '9');
    list->indices[list->num_changes] = strtol(endptr, &endptr, 10);
    assert(endptr[0] == ',');
    endptr++;
    assert(endptr[0] >= '0' && endptr[0] <= '9');
    list->values[list->num_changes] = strtol(endptr, &endptr, 10);
    list->num_changes++;
  }

  assert(*endptr == '\0');
}

void delta_arr_from_uint_arr(const uint32_t *arr, uint32_t num,
                             delta_array_t *list)
{
  list->len = num;
  list->num_changes = 0;
  if(num == 0) return;

  list->indices[0] = 0;
  list->values[0] = arr[0];
  list->num_changes = 1;

  uint32_t i;
  for(i = 1; i < num; i++)
  {
    if(arr[i] != arr[i-1])
    {
      if(list->num_changes == list->changes_capacity)
      {
        list->changes_capacity *= 2;
        size_t mem_uint = list->changes_capacity * sizeof(uint32_t);
        list->indices = realloc(list->indices, mem_uint);
        list->values = realloc(list->values, mem_uint);
      }
      list->indices[list->num_changes] = i;
      list->values[list->num_changes] = arr[i];
      list->num_changes++;
    }
  }
}
