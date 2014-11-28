#ifndef CHROM_POS_LIST_H_
#define CHROM_POS_LIST_H_

// ChromPosOffset coords are 1-based
typedef struct
{
  char *chrom;
  size_t start, end, offset; // always true: start <= end
  bool fw_strand;
} ChromPosOffset;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(chrompos_buf, ChromPosBuffer, ChromPosOffset);

// Sort by length, chrom, strand (fw,rv), start
int chrom_pos_cmp_len(const void *aa, const void *bb);

// Validate a chrom position object
void chrom_pos_validate(const ChromPosOffset *pos);
#define chrom_pos_len(pos) ((pos)->end - (pos)->start)

// Get largest match
// return 1 if there is a unique longest and copies it to pos, otherwise 0
int chrom_pos_list_get_largest(const ChromPosBuffer *buf, ChromPosOffset *pos);

// Parse a string in the form: chr:start-end:strand:offset[,...]
// Return 0 on success, -1 on error
int chrom_pos_list_parse(char *str, ChromPosBuffer *buf);

void chrom_pos_list_sort(ChromPosBuffer *buf);

#endif /* CHROM_POS_LIST_H_ */
