#ifndef FILE_FILTER_H_
#define FILE_FILTER_H_

#include <stdio.h>
#include "string_buffer/string_buffer.h"

// A filter for reading cortex files

typedef struct {
  uint32_t from, into;
} Filter;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(filter_buf, FilterBuffer, Filter);

typedef struct
{
  // If created by opening a file, details are stored in the following strings:
  // `input` => "2:in.ctx:0,3",
  // `path`  => "in.ctx"
  StrBuf input, path;
  uint32_t srcncols; // number of colours in source file/graph
  FilterBuffer filter;
  uint32_t fromncols, intoncols; // max(filter[...].from), max(filter[..].into)
} FileFilter;

// Fetch strings
#define file_filter_input(fltr) ((const char*)(fltr)->input.b)
#define file_filter_path(fltr) ((const char*)(fltr)->path.b)

// Get the number of Filters within a FileFilter
#define file_filter_num(fltr) ((fltr)->filter.len)

// Parse path and create FileFilter, calls die() with msg on error
// call file_filter_set_cols(...) before using this object
void file_filter_open(FileFilter *fltr, const char *path);
// Create a direct filter, ready to be used on return
void file_filter_create_direct(FileFilter *fltr, size_t srcncols, size_t ncols);
// Free memory
void file_filter_close(FileFilter *file);

// If opened with file_filter_open(...), the FileFilter is not ready to be used
// until you have called file_filter_set_cols(...).
// If there is no 'into' filter (e.g. 0,1:in.ctx:4,6 => 0,1 is 'into filter'),
// by default loads into `into_offset..into_offset+N-1`
void file_filter_set_cols(FileFilter *fltr, size_t srcncols, size_t into_offset);

// @param add  Amount to add to each value of intocols
void file_filter_shift_cols(FileFilter *fltr, size_t add);

#define file_filter_intocol(fltr,idx) ((fltr)->filter.b[idx].into)
#define file_filter_fromcol(fltr,idx) ((fltr)->filter.b[idx].from)
#define file_filter_isstdin(fltr) (strcmp((fltr)->path.b,"-") == 0)

#define file_filter_from_ncols(fltr) ((fltr)->fromncols)
#define file_filter_into_ncols(fltr) ((fltr)->intoncols)

uint32_t file_filter_from_ncols_calc(const FileFilter *fltr);
uint32_t file_filter_into_ncols_calc(const FileFilter *fltr);

// Recalculate fromncols, intoncols
// Needs to be called after each modification to the filter
void file_filter_update(FileFilter *fltr);

// Returns true if each colour is loaded directly into the same colour
// 0->0, ... N-1->N-1 for some N
bool file_filter_direct(const FileFilter *fltr);
// Set N to srcncols, dstncols or both
#define file_filter_from_direct(fltr) (file_filter_direct(fltr) && (fltr)->fromncols == (fltr)->srcncols)
#define file_filter_into_direct(fltr,dstncols) (file_filter_direct(fltr) && (fltr)->intoncols == (dstncols))
#define file_filter_full_direct(fltr,dstncols) (file_filter_direct(fltr) && (fltr)->fromncols == (fltr)->srcncols && (fltr)->intoncols == (dstncols))

// Returns true if the specifed filter `fltr` updates colour `col`
bool file_filter_iscolloaded(const FileFilter *fltr, size_t col);

// Print object
void file_filter_status(const FileFilter *fltr, bool writing);

// Clone struct, need to call file_filter_close() to release memory
FileFilter* file_filter_copy(FileFilter *newfltr, const FileFilter *fltr);

/**
 * Set FileFilter to load all samples into a given single colour
 * @param intocols value to set all intocols to
 */
void file_filter_flatten(FileFilter *fltr, size_t intocol);

// Updates filter a using filter b (push a through b: a->b)
// Note: sorts filters in @b
void file_filter_merge(FileFilter *a, FileFilter *b);

// Max from, returns 0 if empty
static inline uint32_t file_filter_max_from(const Filter *fltrs, size_t n) {
  size_t i; uint32_t maxfrom = 0;
  for(i = 0; i < n; i++) maxfrom = MAX2(fltrs[i].from, maxfrom);
  return maxfrom;
}

// Max into, returns 0 if empty
static inline uint32_t file_filter_max_into(const Filter *fltrs, size_t n) {
  size_t i; uint32_t maxinto = 0;
  for(i = 0; i < n; i++) maxinto = MAX2(fltrs[i].into, maxinto);
  return maxinto;
}

static inline int _into_cmp(const void *a, const void *b) {
  int c = cmp(((const Filter*)a)->into, ((const Filter*)b)->into);
  return c ? c : cmp(((const Filter*)a)->from, ((const Filter*)b)->from);
}

static inline int _from_cmp(const void *a, const void *b) {
  int c = cmp(((const Filter*)a)->from, ((const Filter*)b)->from);
  return c ? c : cmp(((const Filter*)a)->into, ((const Filter*)b)->into);
}

#define filters_sort_by_into(fltrs,nfltrs) qsort(fltrs, nfltrs, sizeof(fltrs[0]), _into_cmp)
#define filters_sort_by_from(fltrs,nfltrs) qsort(fltrs, nfltrs, sizeof(fltrs[0]), _from_cmp)


#endif /* FILE_FILTER_H_ */
