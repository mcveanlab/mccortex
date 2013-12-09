#ifndef FILE_FILTER_H_
#define FILE_FILTER_H_

#include "graph_typedef.h"
#include "string_buffer.h"
#include <stdio.h>

// A filter for reading cortex files

typedef struct
{
  // if `orig` was "2:in.ctx:0,3",
  // path would be "in.ctx"
  StrBuf orig, path;
  FILE *fh;
  size_t intocol, ncols, *cols, ncolscap, filencols;
  boolean flatten; // Merge all colours into intocol
  off_t file_size;
  boolean nofilter;
} FileFilter;

#define INIT_FILE_FILTER_MACRO {                                              \
  .orig = {.buff = NULL}, .path = {.buff = NULL}, .fh = NULL,                 \
  .intocol = 0, .ncols = 0, .cols = NULL, .ncolscap = 0, .flatten = false,    \
  .file_size = 0, .nofilter = false}

#define file_filter_outncols(fltr) ((fltr)->flatten ? 1 : (fltr)->ncols)
#define file_filter_intocol(fltr,col) ((fltr)->intocol + (!(fltr)->flatten)*(col))
#define file_filter_fromcol(fltr,col) ((fltr)->cols[i])
#define file_filter_usedcols(fltr) ((fltr)->intocol + file_filter_outncols(fltr))

// Does not read any bytes from file, but does open it
// returns true on success
// on failure will call die (if fatal == true) or return 0 (if fatal == false) 
boolean file_filter_alloc(FileFilter *file, char *path,
                          const char *mode, boolean fatal);

void file_filter_set_cols(FileFilter *fltr, size_t filencols);
void file_filter_update_intocol(FileFilter *fltr, size_t intocol);

void file_filter_close(FileFilter *file);

void file_filter_dealloc(FileFilter *file);

// Print object
void file_filter_status(const FileFilter *fltr);

void file_filter_check_kmer_size(size_t kmer_size, const char *path);

#endif /* FILE_FILTER_H_ */
