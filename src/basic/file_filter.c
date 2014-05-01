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
    file->cols = ctx_malloc(ncolscap * sizeof(*(file->cols)));
    file->ncolscap = ncolscap;
  }
  else if(file->ncolscap < ncolscap) {
    file->cols = ctx_realloc(file->cols, ncolscap * sizeof(*(file->cols)));
    file->ncolscap = ncolscap;
  }
}

// Does not read any bytes from file, but does open it
// returns true on success
// on failure will call die (if fatal == true) or return 0 (if fatal == false)
bool file_filter_alloc(FileFilter *fltr, char *path,
                          const char *mode, bool fatal)
{
  char *path_start, *path_end, path_lchar;

  // Close file if already open
  file_filter_close(fltr);

  if(fltr->orig_path.buff == NULL) strbuf_alloc(&fltr->orig_path, 1024);
  if(fltr->file_path.buff == NULL) strbuf_alloc(&fltr->file_path, 1024);
  strbuf_set(&fltr->orig_path, path);

  file_filter_deconstruct_path(path, &path_start, &path_end);
  fltr->intocol = (size_t)(path_start == path ? 0 : atoi(path));

  path_lchar = *path_end;
  *path_end = '\0';
  strbuf_set(&fltr->file_path, path_start);
  *path_end = path_lchar;

  // Read from stdin in path is "-"
  if(strcmp(fltr->file_path.buff, "-") == 0) {
    fltr->fh = stdin;
    fltr->file_size = -1;
  }
  else
  {
    if((fltr->file_size = futil_get_file_size(fltr->file_path.buff)) == -1) {
      if(fatal) die("Cannot get file size: %s", fltr->file_path.buff);
      else return 0;
    }
    if((fltr->fh = fopen(fltr->file_path.buff, mode)) == NULL) {
      if(fatal) die("Cannot open file: %s", fltr->file_path.buff);
      else return 0;
    }
  }

  return 1;
}

void file_filter_set_cols(FileFilter *fltr, size_t filencols)
{
  ctx_assert(filencols > 0);
  size_t i;
  char *path_start, *path_end;
  file_filter_deconstruct_path(fltr->orig_path.buff, &path_start, &path_end);

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
  if(fltr->intocol != intocol && fltr->intocol != 0)
    warn("Setting load into colour to %zu (%s)", intocol, fltr->orig_path.buff);

  fltr->intocol = intocol;
  for(i = 0; i < fltr->ncols && fltr->cols[i] == i; i++);
  fltr->nofilter = (i == fltr->filencols && fltr->intocol == 0);
}

// Close file
void file_filter_close(FileFilter *fltr)
{
  if(fltr->fh != NULL) { fclose(fltr->fh); fltr->fh = NULL; }
}

void file_filter_dealloc(FileFilter *fltr)
{
  file_filter_close(fltr);
  if(fltr->cols != NULL) { ctx_free(fltr->cols); fltr->cols = NULL, fltr->ncolscap = 0; }
  if(fltr->orig_path.buff != NULL) { strbuf_dealloc(&fltr->orig_path); }
  if(fltr->file_path.buff != NULL) { strbuf_dealloc(&fltr->file_path); }
}

// Print file filter description
void file_filter_status(const FileFilter *fltr)
{
  size_t i;
  const char *file;

  timestamp();
  file = strcmp(fltr->file_path.buff,"-") == 0 ? "STDIN" : fltr->file_path.buff;
  message(" Loading file %s [%zu colour%s]", file, fltr->filencols,
          fltr->filencols != 1 ? "s" : "");
  if(!fltr->nofilter) {
    message(" with colour filter: %zu", fltr->cols[0]);
    for(i = 1; i < fltr->ncols; i++) message(",%zu", fltr->cols[i]);
  }
  size_t into_ncols = file_filter_outncols(fltr);
  if(into_ncols == 1)
    message(" %sinto colour %zu\n", fltr->flatten ? "all " : "", fltr->intocol);
  else
    message(" into colours %zu-%zu\n", fltr->intocol, fltr->intocol+into_ncols-1);
}
