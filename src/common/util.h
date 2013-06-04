
#ifndef UTIL_H_
#define UTIL_H_

int int_cmp(const void *a, const void *b);
int pointer_address_cmp(const void *a, const void *b);

char parse_entire_int(char *str, int *result);
char parse_entire_uint(char *str, unsigned int *result);
char parse_entire_ulong(char *str, unsigned long *result);

// Load a string of numbers into an array. Separator can be any non-numerical
// character. e.g. '1,2,3' '1 2 3'
// Returns 1 if the number of uint parsed == explen, 0 otherwise
boolean parse_uint_liststr(const char *str, uint32_t *arr, uint32_t explen);

size_t count_char(const char *str, char c);

boolean bases_to_integer(const char *arg, size_t *bases);
boolean mem_to_integer(const char *arg, size_t *bytes);

//
// Maths
//

float log_factorial(int number);
float log_factorial_ll(long long number);
unsigned long calculate_mean_ulong(unsigned long *array, unsigned long len);

#endif /* UTIL_H_ */
