#ifndef FILE_UTIL_H_
#define FILE_UTIL_H_

// needed for mode_t used by futil_mkpath(const char *path, mode_t mode)
// and futil_get_file_size(const char* path)
#include <sys/stat.h>
#include "string_buffer/string_buffer.h"


/**
** Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
** Returns 0 on success, -1 on failure
*/
int futil_mkpath(const char *path, mode_t mode);


bool futil_file_exists(const char *file);

bool futil_is_file_readable(const char *file);
// Creates file if it can write
bool futil_is_file_writable(const char *file);
off_t futil_get_file_size(const char* filepath);

#define futil_outpath_str(path) (strcmp(path,"-") == 0 ? "STDOUT" : (path))
#define futil_inpath_str(path) (strcmp(path,"-") == 0 ? "STDIN" : (path))

// Open output file without error if file already exists
FILE *futil_open_output2(const char *path, const char *mode);
gzFile futil_gzopen_output2(const char *path, const char *mode);

// Open and return output file.
// If "-" return stdout, if cannot open die with error message
FILE* futil_open_output(const char *path);
gzFile futil_gzopen_output(const char *path);

// Open and return input file.
// If "-" return stdin, if cannot open die with error message
FILE* futil_open_input(const char *path);
gzFile futil_gzopen_input(const char *path);

// Open a new output file with unused name
bool futil_generate_filename(const char *base_fmt, StrBuf *str);
void futil_get_strbuf_of_dir_path(const char *path, StrBuf *dir);
char* futil_get_current_dir(char abspath[PATH_MAX+1]);

// Case insensitive comparision of path with given extension
bool futil_path_has_extension(const char *path, const char *ext);

// Usage:
//   FILE **tmp_files = futil_create_tmp_files(num_tmp);
// To clear up:
//   for(i = 0; i < num_tmp; i++) fclose(tmp_files[i]);
//   ctx_free(tmp_files);
FILE** futil_create_tmp_files(size_t num_tmp_files);

// Merge temporary files, closes tmp files
void futil_merge_tmp_files(FILE **tmp_files, size_t num_files, FILE *fout);


// This is the same as futil_safe_fread, except it calls return if not `fatal`
#define SAFE_READ(fh,ptr,size,field,path,fatal) {                              \
  size_t _read = fread(ptr,1,size,fh);                                         \
  if(_read != size) {                                                          \
    if(!fatal) return -1;                                                      \
    die("Couldn't read '%s': expected %zu; recieved: %zu; [file: %s]\n",       \
        field, (size_t)size, _read, path);                                     \
  }                                                                            \
}

#define safe_fread(fh,ptr,size,field,path) \
        futil_safe_fread(fh,ptr,size,field,path,__FILE__,__func__,__LINE__)

// fh is where to read from
// ptr is where to load data to
// size is how many bytes to read
// field is the item being read, path is the file being read
// file, line are the position in the codebase
// Moved this to header to help compiler inline
static inline void futil_safe_fread(FILE *fh, void *ptr, size_t size,
                                    const char *field, const char *path,
                                    const char *file, const char *func, int line)
{
  size_t read_nbytes = fread(ptr, 1, size, fh);
  if(read_nbytes != size)
  {
    dief(file, func, line,
             "Couldn't read '%s': expected %zu; recieved: %zu; [file: %s]\n",
             field, size, read_nbytes, path);
  }
}


#endif /* FILE_UTIL_H_ */
