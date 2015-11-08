#ifndef CALL_FILE_READER_H_
#define CALL_FILE_READER_H_

#include "common_buffers.h"

typedef struct
{
  StrBuf txt;
  CharPtrBuffer lines;
  SizeBuffer linelens;
} CallFileEntry;

void call_file_entry_alloc(CallFileEntry *entry);
void call_file_entry_dealloc(CallFileEntry *entry);

// Returns 1 on success 0 on end of file
// dies with error message on bad file
int call_file_read(gzFile gzin, const char *path, CallFileEntry *entry);

#define call_file_num_lines(entry) ((entry)->lines.len)
#define call_file_get_line(entry,i) ((entry)->lines.b[i])
#define call_file_line_len(entry,i) ((entry)->linelens.b[i])

static inline size_t call_file_max_allele_len(const CallFileEntry *centry)
{
  size_t i, max = 0, nlines = call_file_num_lines(centry);
  for(i = 5; i < nlines; i+=2)
    max = MAX2(max, call_file_line_len(centry, i));
  return max;
}

static inline size_t call_file_min_allele_len(const CallFileEntry *centry)
{
  size_t i, min = SIZE_MAX, nlines = call_file_num_lines(centry);
  for(i = 5; i < nlines; i+=2)
    min = MIN2(min, call_file_line_len(centry, i));
  return min;
}

//
// Generally useful functions
//

// bubble_format ? ">bubble." : ">brkpnt."
/**
 * Parse header line from FASTA to fetch call id.
 * Expect ">bubble.<id>." or ">brkpnt.<id>."
 * @param  hdrline  String to read, expect: <leadstr><callid>.
 * @param  leadstr  Expected start of hdrline
 * @return callid or -1 on error
 */
int64_t call_file_get_call_id(const char *hdrline, const char *leadstr);

#endif /* CALL_FILE_READER_H_ */
