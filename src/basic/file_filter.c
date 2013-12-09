#include "global.h"
#include "file_filter.h"
#include "range.h"
#include "file_util.h"

// Get pointers to start and end of actual path
// (\d+:)?path.ctx(:\d+(-\d+)?(,\d+(-\d+)?)*)?
static inline void file_filter_deconstruct_path(char *path,
                                                char **start, char **end)
{
  char *ptr, *c;
  *start = path;
  for(ptr = path; *ptr >= '0' && *ptr <= '9'; ptr++);
  if(ptr > path && *ptr == ':') { ptr++; *start = ptr; }
  // Count backwards to match /:[-,0123456789]*$/
  c = *end = path + strlen(path);
  while(c > (*start)+1) {
    c--;
    if(*c == ':') { *end = c; break; }
    else if(!(*c == ',' || *c == '-' || (*c >= '0' && *c <= '9'))) break;
  }
}

static inline void file_filter_capacity(FileFilter *file, size_t ncolscap)
{
  if(ncolscap == 0) return;
  else if(file->ncolscap == 0) {
    file->cols = malloc2(ncolscap * sizeof(*(file->cols)));
    file->ncolscap = ncolscap;
  }
  else if(file->ncolscap < ncolscap) {
    file->cols = realloc2(file->cols, ncolscap * sizeof(*(file->cols)));
    file->ncolscap = ncolscap;
  }
}

// Does not read any bytes from file, but does open it
// returns true on success
// on failure will call die (if fatal == true) or return 0 (if fatal == false) 
boolean file_filter_alloc(FileFilter *fltr, char *path,
                          const char *mode, boolean fatal)
{
  char *path_start, *path_end, path_lchar;

  // Close file if already open
  if(fltr->fh != NULL) file_filter_close(fltr);

  if(fltr->orig.buff == NULL) strbuf_alloc(&fltr->orig, 1024);
  if(fltr->path.buff == NULL) strbuf_alloc(&fltr->path, 1024);
  strbuf_set(&fltr->orig, path);

  file_filter_deconstruct_path(path, &path_start, &path_end);
  fltr->intocol = (path_start == path ? 0 : atoi(path));

  path_lchar = *path_end;
  *path_end = '\0';
  strbuf_set(&fltr->path, path_start);
  *path_end = path_lchar;

  fltr->file_size = get_file_size(fltr->path.buff);
  if(fltr->file_size == -1) die("Cannot get file size: %s", fltr->path.buff);

  if((fltr->fh = fopen(fltr->path.buff, mode)) == NULL) {
    if(fatal) die("Cannot open file: %s", fltr->path.buff);
    else return 0;
  }

  return 1;
}

void file_filter_set_cols(FileFilter *fltr, size_t filencols)
{
  assert(filencols > 0);
  size_t i;
  char *path_start, *path_end;
  file_filter_deconstruct_path(fltr->orig.buff, &path_start, &path_end);

  fltr->filencols = filencols;

  if(*path_end == ':') {
    fltr->ncols = range_get_num(path_end+1, filencols-1);
    file_filter_capacity(fltr, fltr->ncols);
    range_parse_array(path_end+1, fltr->cols, filencols-1);
    for(i = 0; i < filencols && fltr->cols[i] == i; i++);
    fltr->nofilter = (i == filencols && fltr->intocol == 0);
  }
  else {
    fltr->ncols = filencols;
    file_filter_capacity(fltr, fltr->ncols);
    for(i = 0; i < filencols; i++) fltr->cols[i] = i;
    fltr->nofilter = (fltr->intocol == 0);
  }
}

// Update intocol and nofilter
void file_filter_update_intocol(FileFilter *fltr, size_t intocol)
{
  size_t i;
  fltr->intocol = intocol;
  for(i = 0; i < fltr->ncols && fltr->cols[i] == i; i++);
  fltr->nofilter = (i == fltr->ncols && fltr->intocol == 0);
}

// Close file
void file_filter_close(FileFilter *fltr)
{
  if(fltr->fh != NULL) { fclose(fltr->fh); fltr->fh = NULL; }
}

void file_filter_dealloc(FileFilter *fltr)
{
  file_filter_close(fltr);
  if(fltr->cols != NULL) { free(fltr->cols); fltr->cols = NULL, fltr->ncolscap = 0; }
  if(fltr->orig.buff != NULL) { strbuf_dealloc(&fltr->orig); }
  if(fltr->path.buff != NULL) { strbuf_dealloc(&fltr->path); }
}

// Print file filter description
void file_filter_status(const FileFilter *fltr)
{
  size_t i;
  timestamp(ctx_msg_out);
  message(" Loading file %s", fltr->path.buff);
  if(!fltr->nofilter) {
    message(" with colour filter: %zu", fltr->cols[0]);
    for(i = 1; i < fltr->ncols; i++) message(",%zu", fltr->cols[i]);
  }
  size_t into_ncols = file_filter_outncols(fltr);
  if(into_ncols == 1)
    message(" into colour %zu\n", fltr->intocol);
  else
    message(" into colours %zu-%zu\n", fltr->intocol, fltr->intocol+into_ncols-1);
}

void file_filter_check_kmer_size(size_t kmer_size, const char *path)
{
  const size_t mink = MIN_KMER_SIZE, maxk = MAX_KMER_SIZE;
  if(kmer_size < mink || kmer_size > maxk)
    die("Cannot handle kmer size %zu [%zu-%zu; %s]", kmer_size, mink, maxk, path);
  else if(!(kmer_size & 1))
    die("Kmer size appears to be even! %zu [%s]", kmer_size, path);
}
