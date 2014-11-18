#include "global.h"
#include "call_file_reader.h"
#include "dna.h"

// Print error if an entry has more than 1000 paths
#define MAX_LINES_LIMIT 2004

void call_file_entry_alloc(CallFileEntry *entry)
{
  strbuf_alloc(&entry->txt, 2048);
  char_ptr_buf_alloc(&entry->lines, 64);
}

void call_file_entry_dealloc(CallFileEntry *entry)
{
  strbuf_dealloc(&entry->txt);
  char_ptr_buf_dealloc(&entry->lines);
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
static int fasta_gzread(gzFile gzin, bool ret_on_empty, StrBuf *txt)
{
  size_t linestrt = txt->end;
  // Read past empty lines and comments
  while(strbuf_gzreadline(txt, gzin) > 0 &&
        (txt->b[linestrt] == '#' || txt->b[linestrt] == '\n'))
  {
    // Empty line or comment
    if(ret_on_empty && txt->b[linestrt] == '\n') return 0;
    // Reset line
    txt->b[txt->end = linestrt] = '\0';
  }
  if(txt->end == linestrt) return 0; // no more entries
  if(txt->b[linestrt] != '>') return -1; // bad line
  if(strbuf_gzreadline(txt, gzin) == 0) return -1; // missing sequence
  return 1;
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
  int r = fasta_gzread(gzin, false, txt);
  if(r == 0) return 0;
  else if(r < 0) die("Bad entry [%s]: %s", path, txt->b);

  // Read remainder
  for(i = 0; (r = fasta_gzread(gzin, true, txt)) > 0; i++) {}
  if(r < 0) die("Bad entry [%s]: %s", path, txt->b);

  if(i < 2) die("Too few entries [%s]: %s", path, txt->b);
  if(txt->end < 2 || txt->b[txt->end-1] != '\n' || txt->b[txt->end-2] != '\n')
    die("Calls must end with 2 EOL chars \\n [%s]: %s", path, txt->b);

  // Trim last two end of line characters
  txt->b[txt->end -= 2] = '\0';

  fprintf(stderr, "Got: %s\n", txt->b);

  // Don't do anything that would cause entry->txt to be realloc'd now!
  // We are storing pointers to each line
  char_ptr_buf_reset(&entry->lines);

  ptr = txt->b;
  end = strchr(ptr, '\n');

  while(1)
  {
    char_ptr_buf_add(&entry->lines, ptr);

    if(end) { *end = '\0'; ptr = end+1; end = strchr(ptr, '\n'); }
    else break;

    // Safety measure
    if(entry->lines.len > MAX_LINES_LIMIT)
      die("Too many lines: %zu? [path: %s]: %s", entry->lines.len, path, ptr);
  }

  // Parsing checks
  if(entry->lines.len < 6) die("Fewer than six lines [path: %s]", path);

  // Check lines start with >
  for(i = 0; i < entry->lines.len; i+=2) {
    // printf("a) %s\n", entry->lines.data[i]);
    // printf("b) %s\n", entry->lines.data[i+1]);
    if(entry->lines.data[i][0] != '>')
      die("Line doesn't start with '>' [path: %s]: %s", path, entry->lines.data[i]);
  }

  /*
  // Pull out flanks
  const char *fl5p = entry->lines.data[1];
  size_t fl5plen = call_file_line_len(entry, 1);
  const char *fl3p = entry->lines.data[3];
  size_t fl3plen = call_file_line_len(entry, 3);

  // Reverse complement flank3p
  strbuf_reset(&entry->flank3p_rc);
  strbuf_ensure_capacity(&entry->flank3p_rc, fl3plen);
  dna_revcomp_str(entry->flank3p_rc.b, fl3p, fl3plen);
  entry->flank3p_rc.b[entry->flank3p_rc.len = fl3plen] = '\0';

  // Set shift left / shift right
  size_t num_alleles = entry->lines.len/2 - 2;
  size_t min_left = SIZE_MAX, min_right = SIZE_MAX;
  size_t left = 0, right = 0;

  for(i = 0; i < num_alleles && (min_left || min_right); i++)
  {
    const char *allele = entry->lines.data[i*2+5];
    size_t allelelen = call_file_line_len(entry, i*2+5);
    get_left_right_shift(fl5p, fl5plen, allele, allelelen, fl3p, fl3plen,
                         &left, &right);
    min_left = MIN2(min_left, left);
    min_right = MIN2(min_right, right);
  }

  entry->shift_left = min_left;
  entry->shift_right = min_right;
  */

  return 1;
}
