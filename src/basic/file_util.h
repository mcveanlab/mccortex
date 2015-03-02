#ifndef FILE_UTIL_H_
#define FILE_UTIL_H_

// needed for mode_t used by futil_mkpath(const char *path, mode_t mode)
// and futil_get_file_size(const char* path)
#include <sys/stat.h> // mkdir, mode_t
#include "string_buffer/string_buffer.h"

// Whether to overwrite files without warning (default: false)
bool futil_get_force();
void futil_set_force(bool f);

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

/*!
  Create file and parent directories
  @param path Path to file
  @param mode e.g. O_CREAT | O_EXCL | O_WRONLY | O_APPEND
  @return file descriptor on sucess, -1 if file error and sets errno
 */
int futil_create_file(const char *path, int mode);

/*!
  Create file with write permission, call die() on error.
  If already exists and !futil_get_force() gives error
  @param path does nothing if path is "-" or NULL.
 */
void futil_create_output(const char *path);

/*!
  Open a file and set the buffer to be DEFAULT_IO_BUFSIZE. Call die() if cannot
  open the file.
  @param path If "-" return stdout
  @param mode one of: "r","rw","rw+","a"
 */
FILE* futil_fopen(const char *path, const char *mode);
gzFile futil_gzopen(const char *path, const char *mode);

/*!
  Create file, open and set buffer size. Call die() if cannot create/open.
  If file already exists and !futil_get_force() call die() with error.
  @param path If "-" return stdout
  @param mode one of: "r","rw","rw+","a"
 */
FILE* futil_open_create(const char *path, const char *mode);
gzFile futil_gzopen_create(const char *path, const char *mode);

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
