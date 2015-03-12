#ifndef STR_PARSING_H_
#define STR_PARSING_H_

#include "common_buffers.h"

// Parse a comma separated list e.g. "12,3,12"
// Returns <0 on error, otherwise number of chars used
int comma_list_to_array(const char *str, SizeBuffer *nums);

#endif /* STR_PARSING_H_ */
