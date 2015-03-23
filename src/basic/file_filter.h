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
  // if `orig_path` was "2:in.ctx:0,3",
  // path would be "in.ctx"
  StrBuf input, path;
  uint32_t filencols; // number of colours in file
  FilterBuffer filter;
} FileFilter;

// Fetch strings
#define file_filter_input(fltr) ((const char*)(fltr)->input.b)
#define file_filter_path(fltr) ((const char*)(fltr)->path.b)

// Get the number of Filters within a FileFilter
#define file_filter_num(fltr) ((fltr)->filter.len)

// Parse path and create FileFilter, calls die() with msg on error
void file_filter_open(FileFilter *fltr, const char *path);
// Free memory
void file_filter_close(FileFilter *file);

// FileFilter is not ready to be used until you have called set_cols
// If there is no 'into' filter (e.g. 0,1:in.ctx 0,1 is 'into filter'),
// load into `into_offset..into_offset+N-1`
void file_filter_set_cols(FileFilter *fltr, size_t filencols, size_t into_offset);

// @param add  Amount to add to each value of intocols
void file_filter_shift_cols(FileFilter *fltr, size_t add);

#define file_filter_intocol(fltr,idx) ((fltr)->filter.b[idx].into)
#define file_filter_fromcol(fltr,idx) ((fltr)->filter.b[idx].from)
#define file_filter_isstdin(fltr) (strcmp((fltr)->path.b,"-") == 0)

uint32_t file_filter_from_ncols(const FileFilter *fltr);
uint32_t file_filter_into_ncols(const FileFilter *fltr);

// Returns true if each colour is loaded directly into the same colour
// 0->0, ... N->N where N is filencols
bool file_filter_is_direct(const FileFilter *fltr);

// Returns true if the specifed filter `fltr` updates colour `col`
bool file_filter_iscolloaded(const FileFilter *fltr, size_t col);

// Print object
void file_filter_status(const FileFilter *fltr);

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

#endif /* FILE_FILTER_H_ */
