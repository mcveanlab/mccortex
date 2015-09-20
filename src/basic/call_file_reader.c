#include "global.h"
#include "call_file_reader.h"
#include "dna.h"
#include "util.h"

// Print error if an entry has more than 1000 paths
#define MAX_LINES_LIMIT 2004

void call_file_entry_alloc(CallFileEntry *entry)
{
  strbuf_alloc(&entry->txt, 2048);
  char_ptr_buf_alloc(&entry->lines, 64);
  size_buf_alloc(&entry->linelens, 64);
}

void call_file_entry_dealloc(CallFileEntry *entry)
{
  strbuf_dealloc(&entry->txt);
  char_ptr_buf_dealloc(&entry->lines);
  size_buf_dealloc(&entry->linelens);
}

void call_file_entry_reset(CallFileEntry *entry)
{
  strbuf_reset(&entry->txt);
  char_ptr_buf_reset(&entry->lines);
}

// Read a single line fasta entry
// of the form: ">name\nseq" with optional new line on the end
// returns:
// * 1 on success
// * 0 if no more entries
// * -1 on error
static int fasta_gzread(gzFile gzin, bool ret_on_empty, StrBuf *txt,
                        const char *path)
{
  size_t linestrt = txt->end;
  // Read past empty lines and comments
  while(util_gzcheck(strbuf_gzreadline(txt, gzin), gzin, path) > 0 &&
        (txt->b[linestrt] == '#' || txt->b[linestrt] == '\n'))
  {
    // Empty line or comment
    if(ret_on_empty && txt->b[linestrt] == '\n') return 0;
    // Reset line
    txt->b[txt->end = linestrt] = '\0';
  }
  if(txt->end == linestrt) return 0; // no more entries
  if(txt->b[linestrt] != '>') return -1; // bad line
  size_t s = util_gzcheck(strbuf_gzreadline(txt, gzin), gzin, path);
  return !s ? -1 : 1; // missing sequence if s == 0
}

// Returns 1 on success 0 on end of file
// dies with error message on bad file
int call_file_read(gzFile gzin, const char *path, CallFileEntry *entry)
{
  char *ptr, *end;
  size_t i;
  call_file_entry_reset(entry);
  StrBuf *txt = &entry->txt;

  // Lines starting # are comment lines and should be ignored
  int r = fasta_gzread(gzin, false, txt, path);
  if(r < 0) die("Bad entry [%s]: %s", path, txt->b);
  else if(r == 0) return 0;

  // Read remainder
  for(i = 0; (r = fasta_gzread(gzin, true, txt, path)) > 0; i++) {}
  if(r < 0) die("Bad entry [%s]: %s", path, txt->b);

  if(i < 2) die("Too few entries [%s]: %s", path, txt->b);
  if(txt->end < 2 || txt->b[txt->end-1] != '\n' || txt->b[txt->end-2] != '\n')
    die("Calls must end with 2 EOL chars \\n [%s]: %s", path, txt->b);

  // Trim last two end of line characters
  txt->b[txt->end -= 2] = '\0';

  // printf("W have: '%s'\n", txt->b);

  // Don't do anything that would cause entry->txt to be realloc'd now!
  // We are storing pointers to each line
  char_ptr_buf_reset(&entry->lines);
  size_buf_reset(&entry->linelens);

  ptr = txt->b;

  while(1)
  {
    end = strchr(ptr, '\n');
    size_t len = end ? end - ptr : txt->b + txt->end - ptr;

    char_ptr_buf_push(&entry->lines, &ptr, 1);
    size_buf_push(&entry->linelens, &len, 1);

    if(end) { *end = '\0'; ptr = end+1; }
    else break;

    // Safety measure
    if(entry->lines.len > MAX_LINES_LIMIT)
      die("Too many lines: %zu? [path: %s]: %s", entry->lines.len, path, ptr);
  }

  // Parsing checks
  if(entry->lines.len < 6) die("Fewer than six lines [path: %s]", path);

  // Check lines start with >
  for(i = 0; i < entry->lines.len; i+=2) {
    if(entry->lines.b[i][0] != '>')
      die("Line doesn't start with '>' [path: %s]: %s", path, entry->lines.b[i]);
  }

  return 1;
}
