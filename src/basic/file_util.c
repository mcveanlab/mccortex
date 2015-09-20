#include "global.h"
#include "file_util.h"

#include <libgen.h> // dirname
#include <fcntl.h> // open

bool force_file_overwrite = false;

bool futil_get_force() { return force_file_overwrite; }
void futil_set_force(bool f) { force_file_overwrite = f; }

// Append and remove a byte to the end of the file to force the file
// modification timestamp to update properly
void futil_update_timestamp(const char *path)
{
  struct stat st;
  FILE *fh;

  if(stat(path, &st) != 0 || (fh = fopen(path, "a")) == NULL) {
    warn("Cannot force update file timestamp by appending: %s", path);
  } else {
    fputc('.', fh);
    fclose(fh);
    if(truncate(path, st.st_size) != 0)
      die("Failed to truncate file after adding byte: %s", path);
  }
}


//
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
//

// Returns 0 on success, -1 on failure
static int do_mkdir(const char *path, mode_t mode)
{
  struct stat st;

  // mkdir returns zero on success
  //               nonzero on error, setting errno to EEXIST if dir already exists
  // stat returns nonzero if cannot stat file
  if((mkdir(path, mode) != 0 && errno != EEXIST) ||
     (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)))
    return -1;

  return 0;
}

/**
** Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
** Returns 0 on success, -1 on failure
*/
int futil_mkpath(const char *path, mode_t mode)
{
  char *pp;
  char *sp;
  int status;
  char *copypath = strdup(path);

  status = 0;
  pp = copypath;
  while (status == 0 && (sp = strchr(pp, '/')) != 0)
  {
    if (sp != pp)
    {
      /* Neither root nor double slash in path */
      *sp = '\0';
      status = do_mkdir(copypath, mode);
      *sp = '/';
    }
    pp = sp + 1;
  }
  if(status == 0)
    status = do_mkdir(path, mode);
  free(copypath);
  return status;
}

//
//
//

bool futil_file_exists(const char *file)
{
  return (access(file, F_OK) != -1);
}

bool futil_is_file_readable(const char *file)
{
  FILE *fp = fopen(file, "r");
  if(fp == NULL) return false;
  else fclose(fp);
  return true;
}

// Creates file if it can write
bool futil_is_file_writable(const char *file)
{
  FILE *fp = fopen(file, "a");
  if(fp == NULL) return false;
  else fclose(fp);
  return true;
}

// Returns -1 on failure
off_t futil_get_file_size(const char* filepath)
{
  struct stat st;

  if(strcmp(filepath,"-") == 0) return -1;

  if(stat(filepath, &st) == 0)
    return st.st_size;

  warn("Cannot determine file size: %s [%s]\n", filepath, strerror(errno));

  return -1;
}

/*!
  Create file and parent directories
  @param path Path to file
  @param mode e.g. O_CREAT | O_EXCL | O_WRONLY | O_APPEND
  @param fh   If not NULL, open file handle and return in fh
  @return file descriptor on sucess, -1 if file error and sets errno
 */
int futil_create_file(const char *path, int mode)
{
  int ret, fd;

  ctx_assert(mode & O_CREAT);

  // dirname may modify string, so make copy
  char *pathcpy = strdup(path);
  char *dir = dirname(pathcpy);
  ret = futil_mkpath(dir, 0777);
  free(pathcpy);

  if(ret == -1) return -1;

  if(futil_get_force()) mode &= ~O_EXCL;

  fd = open(path, mode, 0666);

  return fd >= 0 ? fd : -1;
}

/*!
  Check we can create/write to file. If not call die() on error
  @param path does nothing if path is "-" or NULL.
 */
void futil_create_output(const char *path) {
  if(path != NULL && strcmp(path,"-") != 0) {
    int fd = futil_create_file(path, O_CREAT | O_EXCL | O_WRONLY | O_APPEND);
    if(fd < 0) {
      if(errno == EEXIST) die("File already exists: %s", path);
      else die("Cannot write to file: %s [%s]", path, strerror(errno));
    } else {
      close(fd);
    }
  }
}

/*!
  Open a file and set the buffer to be DEFAULT_IO_BUFSIZE. Call die() if cannot
  open the file.
  @param path If "-" return stdout
  @param mode one of: "r","rw","rw+","a"
 */
FILE *futil_fopen(const char *path, const char *mode)
{
  FILE *fout;

  if(path == NULL || strcmp(path,"-") == 0) {
    if(!strcmp(mode,"w")) fout = stdout;
    else if(!strcmp(mode,"r")) fout = stdin;
    else die("Cannot open pipe with mode: %s", mode);
  }
  else if((fout = fopen(path, mode)) == NULL) {
    die("Cannot open file: %s [%s]", futil_outpath_str(path), strerror(errno));
  }

  // Set buffer size
  setvbuf(fout, NULL, _IOFBF, DEFAULT_IO_BUFSIZE);

  return fout;
}

/*!
  @see futil_open()
 */
gzFile futil_gzopen(const char *path, const char *mode)
{
  ctx_assert(strcmp(path, "-") != 0 || strcmp(mode,"w") == 0);
  gzFile gzout = strcmp(path, "-") == 0 ? gzdopen(fileno(stdout), mode)
                                        : gzopen(path, mode);
  if(gzout == NULL)
    die("Cannot open gzfile: %s [%s]", futil_outpath_str(path), strerror(errno));

  // Set buffer size
  #if ZLIB_VERNUM >= 0x1240
    gzbuffer(gzout, DEFAULT_IO_BUFSIZE);
  #endif

  return gzout;
}

/*!
  @abstract Create, open and return output file. Call die() on error.
  @description If file already exists and !futil_get_force() die() with error.
  @param path If "-" return stdout
  @param mode one of: "r","rw","rw+","a"
 */
FILE* futil_fopen_create(const char *path, const char *mode)
{
  ctx_assert(path != NULL);
  futil_create_output(path);
  return futil_fopen(path, mode);
}

/*!
  @see futil_fopen_create()
 */
gzFile futil_gzopen_create(const char *path, const char *mode)
{
  ctx_assert(path != NULL);
  futil_create_output(path);
  return futil_gzopen(path, mode);
}

// Check for errors and close
void futil_fclose(FILE *fh)
{
  if(fh == NULL) return;
  if(fh == stdout || fh == stderr) { fflush(fh); return; }
  if(ferror(fh)) warn("File error: %s [%i]", strerror(errno), errno);
  fclose(fh);
}

// Check for errors and close
void gzutil_fclose(gzFile gz)
{
  if(gz == NULL) return;
  int ecode;
  const char *errstr = gzerror(gz, &ecode);
  if(ecode < 0) warn("GZIP File error: %s [%i]", errstr, ecode);
  gzclose(gz);
}


bool futil_generate_filename(const char *base_fmt, StrBuf *str)
{
  int i;

  for(i = 0; i < 10000; i++)
  {
    strbuf_reset(str);
    strbuf_sprintf(str, base_fmt, i);
    struct stat st;
    if(stat(str->b, &st) != 0) return true;
  }

  return false;
}

// Remember to free the result
void futil_get_strbuf_of_dir_path(const char *path, StrBuf *dir)
{
  char *tmp = strdup(path);
  strbuf_set(dir, dirname(tmp));
  strbuf_append_char(dir, '/');
  ctx_free(tmp);
}

char* futil_get_current_dir(char abspath[PATH_MAX+1])
{
  char cwd[PATH_MAX + 1];
  if(getcwd(cwd, PATH_MAX + 1) != NULL)
    return realpath(cwd, abspath);
  else
    return NULL;
}

// Case insensitive comparision of path with given extension
bool futil_path_has_extension(const char *path, const char *ext)
{
  size_t path_len = strlen(path), ext_len = strlen(ext);
  return (ext_len <= path_len && strcasecmp(path+path_len-ext_len, ext) == 0);
}

// Usage:
//     FILE **tmp_files = futil_create_tmp_files(num_tmp);
// to clear up:
//     for(i = 0; i < num_tmp; i++) fclose(tmp_files[i]);
//     ctx_free(tmp_files);
FILE** futil_create_tmp_files(size_t num_tmp_files)
{
  size_t i;
  StrBuf tmppath;
  strbuf_alloc(&tmppath, 1024);
  FILE **tmp_files = ctx_malloc(num_tmp_files * sizeof(FILE*));

  int r = rand() & ((1<<20)-1);

  for(i = 0; i < num_tmp_files; i++)
  {
    strbuf_reset(&tmppath);
    strbuf_sprintf(&tmppath, "/tmp/cortex.tmp.%i.%zu", r, i);
    if((tmp_files[i] = fopen(tmppath.b, "r+")) == NULL) {
      die("Cannot write temporary file: %s [%s]", tmppath.b, strerror(errno));
    }
    unlink(tmppath.b); // Immediately unlink to hide temp file
  }
  strbuf_dealloc(&tmppath);

  return tmp_files;
}

// Merge temporary files, closes tmp files
void futil_merge_tmp_files(FILE **tmp_files, size_t num_files, FILE *fout)
{
  #define TMP_BUF_SIZE (1<<25) /* 32MB */

  char *data = ctx_malloc(TMP_BUF_SIZE);
  size_t i;
  long len;
  FILE *tmp_file;

  for(i = 0; i < num_files; i++)
  {
    tmp_file = tmp_files[i];
    if(fseek(tmp_file, 0L, SEEK_SET) != 0) die("fseek error");

    while((len = fread(data, 1, TMP_BUF_SIZE, tmp_file)) > 0)
      if(fwrite(data, 1, len, fout) != (unsigned)len)
        die("write error [%s]", strerror(errno));

    if(len < 0) warn("fread error: %s", strerror(errno));
    fclose(tmp_file);
  }

  ctx_free(data);

  #undef TMP_BUF_SIZE
}
