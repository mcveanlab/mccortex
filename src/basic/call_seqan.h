#ifndef CALL_SEQAN_H_
#define CALL_SEQAN_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void multiple_seq_align(char **seqs, int num, char *result, size_t *len);

#ifdef __cplusplus
}
#endif

#endif /* CALL_SEQAN_H_ */
