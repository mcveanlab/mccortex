#ifndef UTIL_H_
#define UTIL_H_

#define util_plural_str(n) ((n) == 1 ? "" : "s")

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
int cmp_charptr(const void *a, const void *b);

bool parse_entire_int(const char *str, int *result);
bool parse_entire_uint(const char *str, unsigned int *result);
bool parse_entire_ulong(const char *str, unsigned long *result);
bool parse_entire_double(const char *str, double *result);
bool parse_entire_size(const char *str, size_t *result);

//
// Bits
//
const uint8_t rev_nibble_table[16];
#define rev_nibble_lookup(x) ({ ctx_assert((unsigned)(x) < 16), rev_nibble_table[(unsigned)(x)]; })

extern const uint8_t nibble_popcount_table[16];

#define byte_popcount(x) (nibble_popcount_table[((x) >> 4) & 0xf] + \
                          nibble_popcount_table[(x) & 0xf])

//
// Numbers to string
//

bool bases_to_integer(const char *arg, size_t *bases);
bool mem_to_integer(const char *arg, size_t *bytes);


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
// Hexidecimal
//

// Generate null terminated string of length num-1
char* hex_rand_str(char *str, size_t num);

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


// sorted_arr must be a sorted array
size_t calc_N50(const size_t *sorted_arr, size_t n, size_t total);

//
// Time
//

// output of form: "10 days 23 hours 59 mins 59 secs"
// extreme: "-2147483647 days 23 hours 59 mins 59 secs",
// so str should be at least 42 bytes long
// returns number of bytes written
size_t seconds_to_str(unsigned long seconds, char *str);

//
// Multi-threading
//

// Run function with given arguments in `nthreads` threads
// Blocks until all jobs finished
void util_run_threads(void *args, size_t nel, size_t elsize,
                      size_t nthreads, void (*func)(void*));

// Increment a uint8_t without overflow
static inline void safe_add_uint8(volatile uint8_t *ptr, uint8_t add)
{
  uint8_t v = *ptr, newv, curr;
  do {
    // Compare and swap returns the value of ptr before the operation
    curr = v;
    newv = MIN2((size_t)UINT8_MAX, (size_t)curr + add);
    if(curr == UINT8_MAX) break;
    v = __sync_val_compare_and_swap(ptr, curr, newv);
  }
  while(v != curr);
}

#define safe_incr_uint8(ptr) safe_add_uint8(ptr,1)

#endif /* UTIL_H_ */
