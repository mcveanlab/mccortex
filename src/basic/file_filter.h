#ifndef FILE_FILTER_H_
#define FILE_FILTER_H_

#include <stdio.h>
#include "string_buffer/string_buffer.h"

// A filter for reading cortex files

typedef struct
{
  // if `orig_path` was "2:in.ctx:0,3",
  // path would be "in.ctx"
  StrBuf orig_path, file_path;
  FILE *fh;
  size_t intocol; // first colour loading into
  size_t ncols; // number of colours being read from file
  size_t *cols, ncolscap; // array of colours to load of length ncols
  size_t filencols; // number of colours in file
  bool flatten; // Merge all colours into intocol
  off_t file_size;
  bool nofilter;
} FileFilter;

#define INIT_FILE_FILTER_MACRO {                                               \
  .orig_path = {.buff = NULL}, .file_path = {.buff = NULL}, .fh = NULL,        \
  .intocol = 0, .ncols = 0, .cols = NULL, .ncolscap = 0, .flatten = false,     \
  .file_size = 0, .nofilter = false}

// returns 1 or ncols [fltr->flatten is 0 or 1]
#define file_filter_outncols(fltr) \
        ((size_t)(fltr)->flatten + (!((fltr)->flatten))*(fltr)->ncols)

#define file_filter_intocol(fltr,col) ((fltr)->intocol + (!(fltr)->flatten)*(col))
#define file_filter_fromcol(fltr,col) ((fltr)->cols[col])
#define file_filter_usedcols(fltr) ((fltr)->intocol + file_filter_outncols(fltr))
#define file_filter_iscolloaded(fltr,col) \
        ((fltr)->intocol<=(col) && (col)<file_filter_usedcols(fltr))

#define file_filter_isstdin(fltr) (strcmp((fltr)->file_path.buff,"-") == 0)


// Does not read any bytes from file, but does open it
// returns true on success
// on failure will call die (if fatal == true) or return 0 (if fatal == false)
bool file_filter_alloc(FileFilter *file, char *path,
                          const char *mode, bool fatal);

// Set nummber of colours in the file
void file_filter_set_cols(FileFilter *fltr, size_t filencols);
// Set which colour to load the first colour into
void file_filter_update_intocol(FileFilter *fltr, size_t intocol);

// Attempt to close file (if open)
void file_filter_close(FileFilter *file);

// Attempt to close file (if open) and free memory
void file_filter_dealloc(FileFilter *file);

// Print object
void file_filter_status(const FileFilter *fltr);

#endif /* FILE_FILTER_H_ */
