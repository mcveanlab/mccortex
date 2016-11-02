/*
 stream_buffer.h
 project: string_buffer
 url: https://github.com/noporpoise/StringBuffer
 author: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain
 Sep 2015
*/

#ifndef _STREAM_BUFFER_HEADER
#define _STREAM_BUFFER_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>

/*
   Generic string buffer functions
*/

#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) _rndup2pow64(x)
  static inline size_t _rndup2pow64(unsigned long long x) {
    // long long >=64 bits guaranteed in C99
    --x; x|=x>>1; x|=x>>2; x|=x>>4; x|=x>>8; x|=x>>16; x|=x>>32; ++x;
    return x;
  }
#endif

static inline void cbuf_capacity(char **buf, size_t *sizeptr, size_t len)
{
  len++; // for nul byte
  if(*sizeptr < len) {
    *sizeptr = ROUNDUP2POW(len);
    if((*buf = realloc(*buf, *sizeptr)) == NULL) {
      fprintf(stderr, "[%s:%i] Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

static inline void cbuf_append_char(char **buf, size_t *lenptr, size_t *sizeptr,
                                    char c)
{
  cbuf_capacity(buf, sizeptr, *lenptr+1);
  (*buf)[(*lenptr)++] = c;
  (*buf)[*lenptr] = '\0';
}

static inline void cbuf_append_str(char **buf, size_t *lenptr, size_t *sizeptr,
                                   const char *str, size_t len)
{
  cbuf_capacity(buf, sizeptr, *lenptr+len);
  memcpy(*buf + *lenptr, str, len);
  *lenptr += len;
  (*buf)[*lenptr] = '\0';
}

static inline void cbuf_chomp(char *buf, size_t *lenptr)
{
  while(*lenptr && (buf[(*lenptr)-1] == '\n' || buf[(*lenptr)-1] == '\r'))
    (*lenptr)--;
  buf[*lenptr] = '\0';
}

/*
   Stream buffer
*/

typedef struct
{
  char *b;
  // begin is index of first char (unless begin >= end)
  // end is index of \0
  // size should be >= end+1 to allow for \0
  // (end-begin) is the number of bytes in buffer
  size_t begin, end, size;
} StreamBuffer;


#define strm_buf_init {.b = NULL, .begin = 0, .end = 0, .size = 0}

#define strm_buf_reset(sb) do { (b)->begin = (b)->end = 1; } while(0)

// Returns NULL if out of memory, @b otherwise
static inline StreamBuffer* strm_buf_alloc(StreamBuffer *b, size_t s)
{
  b->size = (s < 16 ? 16 : s) + 1;
  if((b->b = malloc(b->size)) == NULL) return NULL;
  b->begin = b->end = 1;
  b->b[b->end] = b->b[b->size-1] = 0;
  return b;
}

static inline void strm_buf_dealloc(StreamBuffer *b)
{
  free(b->b);
  memset(b, 0, sizeof(StreamBuffer));
}

static inline StreamBuffer* strm_buf_new(size_t s)
{
  StreamBuffer *b = (StreamBuffer*)malloc(sizeof(StreamBuffer));
  if(b == NULL) return NULL;
  else if(strm_buf_alloc(b,s)) return b;
  else { free(b); return NULL; } /* couldn't malloc */
}

static inline void strm_buf_free(StreamBuffer *cbuf)
{
  free(cbuf->b);
  free(cbuf);
}

// len is the number of bytes you want to be able to store
// Adds one to len for NULL terminating byte
static inline void strm_buf_ensure_capacity(StreamBuffer *cbuf, size_t len)
{
  cbuf_capacity(&cbuf->b, &cbuf->size, len);
}


/*
Unbuffered

fgetc(f)
gzgetc(gz)
fungetc(f,c)
gzungetc(gz,c)
gzread2(gz,buf,len)
fread2(f,buf,len)
gzgets2(gz,buf,len)
fgets2(f,buf,len)
gzreadline(gz,out)
freadline(f,out)
*/

#define ferror2(fh) ferror(fh)
static inline int gzerror2(gzFile gz) { int e; gzerror(gz, &e); return (e < 0); }

// Define read for gzFile and FILE (unbuffered)
// Check ferror/gzerror on return for error
#define fread2(fh,buf,len) fread(buf,1,len,fh)

// gzread returns -1 on error, otherwise number of bytes read
// gzread can only read 2^32 bytes at a time, and returns -1 on error
// Check ferror/gzerror on return for error
static inline size_t gzread2(gzFile gz, void *ptr, size_t len)
{
  size_t nread = 0, n;
  int s;
  while(nread < len) {
    n = len - nread;
    if(n > UINT_MAX) n = UINT_MAX;
    s = gzread(gz, (char*)ptr+nread, n);
    if(s <= 0) break;
    nread += s;
  }
  return nread;
}

// TODO: these should support len > 2^32 by looping fgets/gzgets
#define gzgets2(gz,buf,len) gzgets(gz,buf,(int)(len))
#define fgets2(fh,buf,len) fgets(buf,(int)(len),fh)

// fgetc(f), gzgetc(gz) are already good to go
// fungetc(c,f), gzungetc(c,gz) are already good to go

#define fungetc(c,f) ungetc(c,f)

// Define readline for gzFile and FILE (unbuffered)
// Check ferror/gzerror on return for error
#define _func_readline(name,type_t,__gets) \
  static inline size_t name(type_t file, char **buf, size_t *len, size_t *size)\
  {                                                                            \
    if(*len+1 >= *size) cbuf_capacity(buf, size, *len+1);                      \
    /* Don't read more than 2^32 bytes at once (gzgets limit) */               \
    size_t r = *size-*len > UINT_MAX ? UINT_MAX : *size-*len, origlen = *len;  \
    while(__gets(file, *buf+*len, r) != NULL)                                  \
    {                                                                          \
      *len += strlen(*buf+*len);                                               \
      if((*buf)[*len-1] == '\n') return *len-origlen;                          \
      else *buf = realloc(*buf, *size *= 2);                                   \
      r = *size-*len > UINT_MAX ? UINT_MAX : *size-*len;                       \
    }                                                                          \
    return *len-origlen;                                                       \
  }

_func_readline(gzreadline,gzFile,gzgets2)
_func_readline(freadline,FILE*,fgets2)

// Define skipline
#define _func_skipline(fname,ftype,readc) \
  static inline size_t fname(ftype file)                                       \
  {                                                                            \
    int c;                                                                     \
    size_t skipped_bytes = 0;                                                  \
    while((c = readc(file)) != -1) {                                           \
      skipped_bytes++;                                                         \
      if(c == '\n') break;                                                     \
    }                                                                          \
    return skipped_bytes;                                                      \
  }

_func_skipline(gzskipline,gzFile,gzgetc)
_func_skipline(fskipline,FILE*,fgetc)

/* Buffered */

/*
fgetc_buf(f,in)
gzgetc_buf(gz,in)
ungetc_buf(c,in)
fread_buf(f,ptr,len,in)
gzread_buf(f,ptr,len,in)
gzreadline_buf(gz,in,out)
freadline_buf(f,in,out)
*/

// __read is either gzread2 or fread2
// offset of 1 so we can unget at least one char
// Beware: read-in buffer is not null-terminated
// Returns fail on error
#define _READ_BUFFER(file,in,__read) do                                        \
{                                                                              \
  (in)->end = 1+__read(file,(in)->b+1,(in)->size-1);                           \
  (in)->begin = 1;                                                             \
} while(0)

// Define getc for gzFile and FILE (buffered)
// Returns character or -1 at EOF / error
// Check ferror/gzerror on return for error
#define _func_getc_buf(fname,type_t,__read)                                    \
  static inline int fname(type_t file, StreamBuffer *in)                       \
  {                                                                            \
    if(in->begin >= in->end) {                                                 \
      _READ_BUFFER(file,in,__read);                                            \
      return in->begin >= in->end ? -1 : in->b[in->begin++];                   \
    }                                                                          \
    return in->b[in->begin++];                                                 \
  }

_func_getc_buf(fgetc_buf,FILE*,fread2)
_func_getc_buf(gzgetc_buf,gzFile,gzread2)

// Define ungetc for buffers
// returns c if successful, otherwise -1
static inline int ungetc_buf(int c, StreamBuffer *in)
{
  if(in->begin == 0) {
    if(in->end == 0) {
      in->b[0] = (char)c;
      in->end = 1;
      return c;
    }
    else return -1;
  }
  in->b[--(in->begin)] = (char)c;
  return c;
}

#define fungetc_buf(c,in) ungetc_buf(c,in)
#define gzungetc_buf(c,in) ungetc_buf(c,in)

// Returns number of bytes read
// Check ferror/gzerror on return for error
#define _func_read_buf(fname,type_t,__read)                                    \
  static inline size_t fname(type_t file, void *ptr, size_t len, StreamBuffer *in)\
  {                                                                            \
    size_t remaining = len, next;                                              \
    if(len > 2*in->size) { /* first empty buffer then read directly */         \
      next = in->end - in->begin;                                              \
      memcpy(ptr, in->b+in->begin, in->end-in->begin);                         \
      in->begin = in->end; ptr = (char*)ptr + next; remaining -= next;         \
      remaining -= __read(file,ptr,remaining);                                 \
    }                                                                          \
    else {                                                                     \
      while(1) {                                                               \
        if(remaining == 0) { break; }                                          \
        if(in->begin >= in->end) {                                             \
          _READ_BUFFER(file,in,__read);                                        \
          if(in->begin >= in->end) { break; }                                  \
        }                                                                      \
        next = in->end - in->begin;                                            \
        if(remaining < next) next = remaining;                                 \
        memcpy(ptr, in->b+in->begin, next);                                    \
        in->begin += next; ptr = (char*)ptr + next; remaining -= next;         \
      }                                                                        \
    }                                                                          \
    return len - remaining;                                                    \
  }

_func_read_buf(gzread_buf,gzFile,gzread2)
_func_read_buf(fread_buf,FILE*,fread2)

// Define readline for gzFile and FILE (buffered)
// Check ferror/gzerror on return for error
#define _func_readline_buf(fname,type_t,__read)                                \
  static inline size_t fname(type_t file, StreamBuffer *in,                    \
                             char **buf, size_t *len, size_t *size)            \
  {                                                                            \
    if(in->begin >= in->end) { _READ_BUFFER(file,in,__read); }                 \
    size_t offset, buffered, total_read = 0;                                   \
    while(in->end > in->begin)                                                 \
    {                                                                          \
      for(offset = in->begin; offset < in->end && in->b[offset++] != '\n'; ){} \
      buffered = offset - in->begin;                                           \
      cbuf_capacity(buf, size, (*len)+buffered);                               \
      memcpy((*buf)+(*len), in->b+in->begin, buffered);                        \
      *len += buffered;                                                        \
      in->begin = offset;                                                      \
      total_read += buffered;                                                  \
      if((*buf)[*len-1] == '\n') break;                                        \
      _READ_BUFFER(file,in,__read);                                            \
    }                                                                          \
    (*buf)[*len] = 0;                                                          \
    return total_read;                                                         \
  }

_func_readline_buf(gzreadline_buf,gzFile,gzread2)
_func_readline_buf(freadline_buf,FILE*,fread2)

// Define buffered skipline
// Check ferror/gzerror on return for error
#define _func_skipline_buf(fname,ftype,__read)                                 \
  static inline size_t fname(ftype file, StreamBuffer *in)                     \
  {                                                                            \
    if(in->begin >= in->end) { _READ_BUFFER(file,in,__read); }                 \
    size_t offset, skipped_bytes = 0;                                          \
    while(in->end > in->begin)                                                 \
    {                                                                          \
      for(offset = in->begin; offset < in->end && in->b[offset++] != '\n'; )   \
      skipped_bytes += offset - in->begin;                                     \
      in->begin = offset;                                                      \
      if(in->b[offset-1] == '\n') break;                                       \
      _READ_BUFFER(file,in,__read);                                            \
    }                                                                          \
    return skipped_bytes;                                                      \
  }

_func_skipline_buf(gzskipline_buf,gzFile,gzread2)
_func_skipline_buf(fskipline_buf,FILE*,fread2)

// Define buffered gzgets_buf, fgets_buf

// Reads upto len-1 bytes (or to the first \n if first) into str
// Adds null-terminating byte
// returns pointer to `str` or NULL if EOF
// Check ferror/gzerror on return for error
#define _func_gets_buf(fname,type_t,__read)                                    \
  static inline char* fname(type_t file, StreamBuffer *in, char *str, size_t len)\
  {                                                                            \
    if(len == 0) return NULL;                                                  \
    if(len == 1) {str[0] = 0; return str; }                                    \
    if(in->begin >= in->end) { _READ_BUFFER(file,in,__read); }                 \
    size_t i, buffered, limit, total_read = 0, remaining = len-1;              \
    while(in->end > in->begin)                                                 \
    {                                                                          \
      limit = (in->begin+remaining < in->end ? in->begin+remaining : in->end); \
      for(i = in->begin; i < limit && in->b[i++] != '\n'; );                   \
      buffered = i - in->begin;                                                \
      memcpy(str+total_read, in->b+in->begin, buffered);                       \
      in->begin += buffered;                                                   \
      total_read += buffered;                                                  \
      remaining -= buffered;                                                   \
      if(remaining == 0 || str[total_read-1] == '\n') break;                   \
      _READ_BUFFER(file,in,__read);                                            \
    }                                                                          \
    str[total_read] = 0;                                                       \
    return total_read == 0 ? NULL : str;                                       \
  }

_func_gets_buf(gzgets_buf,gzFile,gzread2)
_func_gets_buf(fgets_buf,FILE*,fread2)


// Buffered ftell/gztell, fseek/gzseek

// ftell
static inline off_t ftell_buf(FILE *fh, StreamBuffer *strm) {
  return ftell(fh) - (strm->end - strm->begin);
}

static inline off_t gztell_buf(gzFile gz, StreamBuffer *strm) {
  return gztell(gz) - (strm->end - strm->begin);
}

// fseek
static inline long fseek_buf(FILE *fh, off_t offset, int whence,
                             StreamBuffer *strm)
{
  // fh points to byte after end of buffer
  // ftell(fh) returns <0 on error
  off_t n = strm->end-strm->begin;
  off_t t = ftell(fh), s = t - n;
  if(t >= 0 && whence == SEEK_CUR && offset >= 0 && offset < n) {
    strm->begin += offset;
    return 0;
  } else if(t >= 0 && whence == SEEK_SET && s <= offset && offset < t) {
    strm->begin += offset-s;
    return 0;
  } else {
    strm->begin = strm->end = 1; // wipe buffer
    return fseek(fh, offset, whence);
   }
}

static inline long gzseek_buf(gzFile gz, off_t offset, int whence,
                              StreamBuffer *strm)
{
  // gz points to byte after end of buffer
  // gztell(gz) returns <0 on error
  off_t n = strm->end-strm->begin;
  off_t t = gztell(gz), s = t - n;
  if(t >= 0 && whence == SEEK_CUR && offset >= 0 && offset < n) {
    strm->begin += offset;
    return 0;
  } else if(t >= 0 && whence == SEEK_SET && s <= offset && offset < t) {
    strm->begin += offset-s;
    return 0;
  } else {
    strm->begin = strm->end = 1; // wipe buffer
    return gzseek(gz, offset, whence);
   }
}

/*
 Output (unbuffered)

fputc2(fh,c)
gzputc2(gz,c)
fputs2(fh,c)
gzputs2(gz,c)
fprintf(fh,fmt,...) / gzprintf(gz,fmt,...) already useable
fwrite2(fh,ptr,len)
gzwrite2(gz,ptr,len)
*/

#define fputc2(fh,c) fputc(c,fh)
#define gzputc2(gz,c) gzputc(gz,c)

#define fputs2(fh,str) fputs(str,fh)
#define gzputs2(gz,str) gzputs(gz,str)

#define fwrite2(fh,ptr,len) fwrite(ptr,len,1,fh)
#define gzwrite2(gz,ptr,len) gzwrite(gz,ptr,len)

/*
 Output (buffered)

// TODO
fputc_buf(fh,buf,c)
gzputc_buf(gz,buf,c)
fputs_buf(fh,buf,str)
gzputs_buf(gz,buf,str)
fprintf_buf(fh,buf,fmt,...)
gzprintf_buf(gz,buf,fmt,...)
fwrite_buf(fh,buf,ptr,len)
gzwrite_buf(gz,buf,ptr,len)
strm_buf_flush(fh,buf)
strm_buf_gzflush(gz,buf)
*/


#endif
