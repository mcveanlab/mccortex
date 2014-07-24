#include "global.h"
#include "file_filter.h"
#include "range.h"
#include "util.h"
#include "file_util.h"

void file_filter_ensure_capacity(FileFilter *fltr, size_t size)
{
  if(fltr->capacity < size) {
    fltr->capacity = roundup2pow(size);
    fltr->filter = ctx_reallocarray(fltr->filter, fltr->capacity, sizeof(Filter));
  }
}

// Get pointers to start and end of actual path
// (\d+:)?path.ctx(:\d+(-\d+)?(,\d+(-\d+)?)*)?
static inline void file_filter_deconstruct_path(const char *path,
                                                const char **start,
                                                const char **end)
{
  const char *ptr, *c;
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

// Parse path and create FileFilter, calls die() with msg on error
// fltr should be zero'd before call
void file_filter_open(FileFilter *fltr, const char *path)
{
  const char *path_start, *path_end;
  size_t path_len;

  // Duplicate input string and file path
  strbuf_set(&fltr->input, path);
  file_filter_deconstruct_path(path, &path_start, &path_end);
  path_len = path_end - path_start;
  strbuf_ensure_capacity(&fltr->path, path_len);
  memcpy(fltr->path.b, path_start, path_len);
  fltr->path.b[fltr->path.end = path_len] = '\0';
}

// Close file, release memory
void file_filter_close(FileFilter *fltr)
{
  strbuf_dealloc(&fltr->input);
  strbuf_dealloc(&fltr->path);
  ctx_free(fltr->filter);
  memset(fltr, 0, sizeof(FileFilter));
}

// Return max colours we load from, or 0 if no colours
uint32_t file_filter_from_ncols(const FileFilter *fltr)
{
  uint32_t i, max = 0;
  for(i = 0; i < fltr->ncols; i++)
    max = MAX2(max, fltr->filter[i].from+1);
  return max;
}

// Return max colours we load into, or 0 if no colours
uint32_t file_filter_into_ncols(const FileFilter *fltr)
{
  uint32_t i, max = 0;
  for(i = 0; i < fltr->ncols; i++)
    max = MAX2(max, fltr->filter[i].into+1);
  return max;
}

void file_filter_set_cols(FileFilter *fltr, size_t filencols)
{
  size_t i;
  const char *path_start, *path_end;
  file_filter_deconstruct_path(fltr->input.b, &path_start, &path_end);

  fltr->orig_first_col = (size_t)(path_start == fltr->input.b ? 0 : atoi(fltr->input.b));

  fltr->ncols = filencols;
  if(*path_end == ':') {
    int s = range_get_num(path_end+1, filencols-1);
    if(s == -1) die("Invalid filter path: %s", fltr->input.b);
    fltr->ncols = s;
  }

  fltr->filencols = filencols;
  file_filter_ensure_capacity(fltr, fltr->ncols);

  if(*path_end == ':')
  {
    size_t *tmp = ctx_calloc(fltr->ncols, sizeof(size_t));
    if(range_parse_array(path_end+1, tmp, filencols-1) == -1)
      die("Invalid filter path: %s", fltr->input.b);
    for(i = 0; i < fltr->ncols; i++)
      fltr->filter[i].from = tmp[i];
    ctx_free(tmp);
  }
  else {
    for(i = 0; i < fltr->ncols; i++)
      fltr->filter[i].from = i;
  }

  for(i = 0; i < fltr->ncols; i++)
    fltr->filter[i].into = fltr->orig_first_col + i;
}

// @add amount to add to each value of intocols
void file_filter_shift_cols(FileFilter *fltr, size_t add)
{
  size_t i;
  for(i = 0; i < fltr->ncols; i++) fltr->filter[i].into += add;
}

// Returns true if each colour is loaded directly into the same colour
// 0->0, ... N->N where N is filencols
bool file_filter_is_direct(const FileFilter *fltr)
{
  size_t i;
  if(fltr->ncols != fltr->filencols) return false;
  for(i = 0; i < fltr->ncols && fltr->filter[i].from == fltr->filter[i].into; i++);
  return (i == fltr->ncols);
}

// Returns true if the specifed filter `fltr` updates colour `col`
bool file_filter_iscolloaded(const FileFilter *fltr, size_t col)
{
  size_t i;
  for(i = 0; i < fltr->ncols; i++) {
    if(file_filter_intocol(fltr, i) == col)
      return true;
  }
  return false;
}

// Print file filter description
void file_filter_status(const FileFilter *fltr)
{
  size_t i;

  pthread_mutex_lock(&ctx_biglock);
  timestamp();
  message("[FileFilter] Loading file %s [%u colour%s]", file_filter_path(fltr),
          fltr->filencols, util_plural_str(fltr->filencols));

  if(!file_filter_is_direct(fltr)) {
    message(" with filter: %u->%u", fltr->filter[0].from, fltr->filter[0].into);
    for(i = 1; i < fltr->ncols; i++)
      message(",%u->%u", fltr->filter[i].from, fltr->filter[i].into);
  }
  message("\n");
  pthread_mutex_unlock(&ctx_biglock);
}

// Copy FileFilter src into dst
FileFilter* file_filter_copy(FileFilter *dst, const FileFilter *src)
{
  strbuf_set_buff(&dst->input, &src->input);
  strbuf_set_buff(&dst->path, &src->path);
  file_filter_ensure_capacity(dst, src->ncols);
  memcpy(dst->filter, src->filter, src->ncols * sizeof(Filter));
  dst->ncols = src->ncols;
  dst->filencols = src->filencols;
  return dst;
}

// @intocols value to set all intocols to
void file_filter_flatten(FileFilter *fltr, size_t intocol)
{
  size_t i;
  ctx_assert(fltr->filter != NULL);
  for(i = 0; i < fltr->ncols; i++)
    fltr->filter[i].into = intocol;
}

void file_filter_add(FileFilter *fltr, uint32_t from, uint32_t into)
{
  file_filter_ensure_capacity(fltr, fltr->ncols+1);
  fltr->filter[fltr->ncols++] = (Filter){.from = from, .into = into};
}

//
// Sort array of Filter structs
//
static inline int _into_cmp(const void *a, const void *b) {
  return (long)((const Filter*)a)->into - ((const Filter*)b)->into;
}

static inline int _from_cmp(const void *a, const void *b) {
  return (long)((const Filter*)a)->from - ((const Filter*)b)->from;
}

#define filters_sort_by_into(fltr) qsort(fltr->filter, fltr->ncols, sizeof(fltr->filter[0]), _into_cmp)
#define filters_sort_by_from(fltr) qsort(fltr->filter, fltr->ncols, sizeof(fltr->filter[0]), _from_cmp)

// Updates filter a using filter b (push a through b: a->b)
// Note: sorts filters in @b
void file_filter_merge(FileFilter *a, FileFilter *b)
{
  // Update a->intocols[i] = b->intocols[j] where a->intocols[i]==b->fromcols[j]
  // if no value of j such that a->intocols[i]==b->fromcols[j], remove j

  filters_sort_by_into(a);
  filters_sort_by_from(b);

  size_t i, j, k, n = a->ncols, m = b->ncols;

  for(i = j = 0; i < n; i++) {
    while(j < m && a->filter[i].into > b->filter[j].from) j++;
    for(k = j; k < m && a->filter[i].into == b->filter[k].from; k++)
      file_filter_add(a, a->filter[i].from, b->filter[k].into);
  }

  memmove(a->filter, a->filter+n, a->ncols - n);
  a->ncols -= n;
}
