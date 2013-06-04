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

boolean test_file_readable(const char *file);
// Creates file if it can write
boolean test_file_writable(const char *file);
off_t get_file_size(const char* filepath);

// Open a new output file with unused name
StrBuf* file_reader_generate_filename(const char *base_fmt);
StrBuf *file_reader_get_strbuf_of_dir_path(char *path);
char* file_reader_get_current_dir(char abspath[PATH_MAX+1]);

#endif /* FILE_UTIL_H_ */
