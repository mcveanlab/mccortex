#ifndef RANGE_H_
#define RANGE_H_

uint32_t range_get_num(const char *str, uint32_t range_max);
void range_parse_array(const char *str, uint32_t *arr, uint32_t range_max);

#endif /* RANGE_H_ */
