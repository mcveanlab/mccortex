#ifndef RANGE_H_
#define RANGE_H_

/*
 * Valid ranges are:
 *   *
 *   1
 *   3
 *   2-4
 *   1,1-3,2
 */

/**
 * Parse range string and return number of items
 *
 * @return number of items in range, or -1 if there is a syntax error
 */
int range_get_num(const char *str, size_t range_max);

/**
 * Parse range string into array arr
 *
 * @param str nul terminated string to parse
 * @param arr place parsed array here
 * @param range_max max value permitted in the array
 * @return 0 on success, -1 on error
 */
int range_parse_array(const char *str, size_t *arr, size_t range_max);

/**
 * Parse range into array arr, filling array to ensure exactly a given number
 * of entries. If empty, array is filled 0..num_entries-1, if only one entry,
 * array is filled with same entry num_entries times
 *
 * @param str nul terminated string to parse
 * @param arr place parsed array here
 * @param range_max max value permitted in the array
 * @param num_entries Force exactly `num_entries` to be placed in `arr`
 * @return 0 on success, -1 on error
 */
int range_parse_array_fill(const char *str, size_t *arr,
                           size_t range_max, size_t num_entries);

#endif /* RANGE_H_ */
