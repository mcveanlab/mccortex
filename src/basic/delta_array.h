
#ifndef DELTA_ARR_H_
#define DELTA_ARR_H_

#include <zlib.h>

typedef struct
{
  uint32_t *indices, *values, len;
  size_t num_changes, changes_capacity;
  // unpacked into
  uint32_t *arr, arr_capacity;
} DeltaArray;

// 100,6,20,7,52,8
//  100,6 => 100 long, arr[0] = 6
//  20,7  => arr[20] = 7
//  52,8  => arr[52] = 8

// 0,0 for list of length zero

void delta_arr_alloc(DeltaArray *list);
void delta_arr_dealloc(DeltaArray *list);

void delta_arr_reverse(DeltaArray *list);
void delta_array_unpack(DeltaArray *list);

void delta_arr_from_str(const char *str, DeltaArray *list);
void delta_arr_from_uint_arr(const uint32_t *arr, uint32_t num,
                             DeltaArray *list);

void delta_arr_print(DeltaArray *list, FILE *fh);
void delta_arr_gzprint(DeltaArray *list, gzFile gz);

#endif /* DELTA_ARR_H_ */
