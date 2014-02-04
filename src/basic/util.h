#ifndef UTIL_H_
#define UTIL_H_

// comparison returns:
//   negative iff a < b
//          0 iff a == b
//   positive iff a > b

int cmp_int(const void *a, const void *b);
int cmp_long(const void *a, const void *b);
int cmp_float(const void *a, const void *b);
int cmp_double(const void *a, const void *b);
int cmp_uint32(const void *a, const void *b);
int cmp_uint64(const void *a, const void *b);
int cmp_size(const void *a, const void *b);
int cmp_ptr(const void *a, const void *b);

char parse_entire_int(const char *str, int *result);
char parse_entire_uint(const char *str, unsigned int *result);
char parse_entire_ulong(const char *str, unsigned long *result);
char parse_entire_double(const char *str, double *result);
char parse_entire_size(const char *str, size_t *result);

size_t count_char(const char *str, char c);

boolean bases_to_integer(const char *arg, size_t *bases);
boolean mem_to_integer(const char *arg, size_t *bytes);


unsigned int num_of_digits(unsigned long num);

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27 bytes
// returns pointer to result
char* ulong_to_str(unsigned long num, char* result);

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('-9,223,372,036,854,775,808')+1 = 27 bytes
char* long_to_str(long num, char* result);

// result must be long enough for result + 1 ('\0').
// Max length required is: 26+1+decimals+1 = 28+decimals bytes
// strlen('-9,223,372,036,854,775,808') = 27
// strlen('.') = 1
// +1 for \0
char* double_to_str(double num, int decimals, char* str);

// str must be 26 + 3 + 1 + num decimals + 1 = 31+decimals bytes
// breakdown:
// strlen('18,446,744,073,709,551,615') = 26
// strlen(' GB') = 3
// strlen('.') = 1
// +1 for '\0'
char* bytes_to_str(unsigned long num, int decimals, char* str);

// Number to string using G to mean 10^9, M to mean 10^6 etc
char* num_to_str(double num, int decimals, char* str);

//
// Floats and Doubles
//

// http://www.devx.com/tips/Tip/42853
// static inline int my_isnan(double x) {
//   volatile double temp = x;
//   return temp != x;
// }
// static inline int my_isinf(double x) {
//   volatile double temp = x;
//   return ((temp == x) && ((temp - x) != 0.0));
// }

//
// Maths
//

float log_factorial(unsigned int number);
float log_factorial_ll(unsigned long long number);
unsigned long calculate_mean_ulong(unsigned long *array, unsigned long len);

float find_hist_median(const uint64_t *arr, size_t arrlen, size_t sum);

uint32_t calc_GCD(uint32_t a, uint32_t b);

//
// Time
//

// output of form: "10 days 23 hours 59 mins 59 secs"
// extreme: "-2147483647 days 23 hours 59 mins 59 secs",
// so str should be at least 42 bytes long
// returns number of bytes written
size_t seconds_to_str(unsigned long seconds, char *str);

#endif /* UTIL_H_ */
