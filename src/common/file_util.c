
#include "global.h"
#include "file_util.h"

#include <libgen.h> // dirname
#include <errno.h>

// Returns 1 on success, 0 on failure
// Sets errno to ENOTDIR if already exists but is not directory
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
static char ensure_dir_exists(const char *path, mode_t mode)
{
  struct stat st;

  if(stat(path, &st) != 0)
  {
    // Directory does not exist
    return mkdir(path, mode) == 0 ? 1 : 0;
  }
  else if(!S_ISDIR(st.st_mode))
  {
    errno = ENOTDIR;
    return 0;
  }

  return 1;
}

// mkpath - ensure all directories in path exist
// Returns 1 on success, 0 on failure
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
char mkpath(const char *path, mode_t mode)
{
  char *copypath = strdup(path);

  size_t i = 0, j = 0;
  char status = 1;

  while(1)
  {
    while(path[i] == '.' || path[i] == '/') i++;
    j = i;

    while(path[j] != '.' && path[j] != '/' && path[j] != '\0') j++;
    if(i == j) break;

    char tmp = copypath[j];
    copypath[j] = '\0';

    if(!(status = ensure_dir_exists(copypath, mode))) break;
    if(tmp == '\0') break;

    copypath[j] = tmp;
    i = j + 1;
  }

  free(copypath);
  return status;
}

boolean test_file_readable(const char *file)
{
  FILE *fp = fopen(file, "r");
  if(fp == NULL) return false;
  else fclose(fp);
  return true;
}

// Creates file if it can write
boolean test_file_writable(const char *file)
{
  FILE *fp = fopen(file, "w");
  if(fp == NULL) return false;
  else fclose(fp);
  return true;
}

// Returns -1 on failure
off_t get_file_size(const char* filepath)
{
  struct stat st;

  if (stat(filepath, &st) == 0)
      return st.st_size;

  warn("Cannot determine size of %s: %s\n", filepath, strerror(errno));

  return -1;
}

boolean file_reader_generate_filename(const char *base_fmt, StrBuf *str)
{
  int i;

  for(i = 0; i < 10000; i++)
  {
    strbuf_reset(str);
    strbuf_sprintf(str, base_fmt, i);
    struct stat st;
    if(stat(str->buff, &st) != 0) return true;
  }

  return false;
}

// Remember to free the result
void file_reader_get_strbuf_of_dir_path(const char *path, StrBuf *dir)
{
  char *tmp = strdup(path);
  strbuf_set(dir, dirname(tmp));
  strbuf_append_char(dir, '/');
  free(tmp);
}

char* file_reader_get_current_dir(char abspath[PATH_MAX+1])
{
  char cwd[PATH_MAX + 1];
  if(getcwd(cwd, PATH_MAX + 1) != NULL)
    return realpath(cwd, abspath);
  else
    return NULL;
}
