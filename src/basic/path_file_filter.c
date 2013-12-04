#include "global.h"
#include "path_file_filter.h"
#include "path_format.h"

const PathFileHeader INIT_PATH_FILE_HDR = INIT_PATH_FILE_HDR_MACRO;
const PathFileReader INIT_PATH_READER = INIT_PATH_READER_MACRO;

int path_file_open(PathFileReader *file, char *path, boolean fatal)
{
  return path_file_open2(file, path, fatal, "r");
}

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new PathFileReader and returns 1
int path_file_open2(PathFileReader *file, char *path, boolean fatal,
                     const char *mode)
{
  PathFileHeader *hdr = &file->hdr;
  FileFilter *fltr = &file->fltr;

  if(!file_filter_alloc(fltr, path, mode, fatal)) return 0;
  setvbuf(fltr->fh, NULL, _IOFBF, CTX_BUF_SIZE);

  file->hdr_size = paths_file_read_header(fltr->fh, hdr, fatal, fltr->path.buff);
  if(file->hdr_size == -1) return -1;

  file_filter_set_cols(fltr, hdr->num_of_cols-1);

  // File header checks

  return 1;
}

// Close file
void path_file_close(PathFileReader *file)
{
  file_filter_close(&file->fltr);
}

// calls file_filter_dealloc which will close file if needed
void path_file_dealloc(PathFileReader *file)
{
  file_filter_dealloc(&file->fltr);
  paths_header_dealloc(&file->hdr);
}
