#ifndef RANGE_H_
#define RANGE_H_

// Returns number of items in range, or -1 if there is a syntax error
int range_get_num(const char *str, size_t range_max);

// Parse range into array arr
// Returns 0 on success, -1 on error
int range_parse_array(const char *str, size_t *arr, size_t range_max);

#endif /* RANGE_H_ */
