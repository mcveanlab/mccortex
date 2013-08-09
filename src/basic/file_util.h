#ifndef FILE_UTIL_H_
#define FILE_UTIL_H_

// needed for mode_t used by mkpath(const char *path, mode_t mode)
// and get_file_size(const char* path)
#include <sys/stat.h>
#include "string_buffer.h"

// mkpath - ensure all directories in path exist
// Returns 1 on success, 0 on failure
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
char mkpath(const char *path, mode_t mode);

boolean file_exists(const char *file);

boolean test_file_readable(const char *file);
// Creates file if it can write
boolean test_file_writable(const char *file);
off_t get_file_size(const char* filepath);

// Open a new output file with unused name
boolean futil_generate_filename(const char *base_fmt, StrBuf *str);
void futil_get_strbuf_of_dir_path(const char *path, StrBuf *dir);
char* futil_get_current_dir(char abspath[PATH_MAX+1]);

#define safe_fread(fh,ptr,size,field,path) \
        futil_safe_fread(fh,ptr,size,field,path,__FILE__,__LINE__)

void futil_safe_fread(FILE *fh, void *ptr, size_t size,
                      const char* field, const char *path,
                      char *file, int line);

size_t stream_skip(FILE *fh, size_t skip);

#define SAFE_READ(fh,ptr,size,field,path,fatal) {                          \
  size_t _read = fread(ptr,1,size,fh);                                     \
  if(_read != size) {                                                      \
    if(!fatal) return -1;                                                  \
    die("Couldn't read '%s': expected %zu; recieved: %zu; [file: %s]\n",   \
        field, (size_t)size, _read, path);                                 \
  }                                                                        \
}

#endif /* FILE_UTIL_H_ */
