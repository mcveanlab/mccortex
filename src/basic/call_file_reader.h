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
#define call_file_get_line(entry,i) ((entry)->lines.data[i])
#define call_file_line_len(entry,i) ((entry)->linelens.data[i])

#endif /* CALL_FILE_READER_H_ */
