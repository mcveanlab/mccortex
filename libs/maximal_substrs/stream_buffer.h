/*
 stream_buffer.c
 project: string_buffer
 url: https://github.com/noporpoise/StringBuffer
 author: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain
 Jan 2014
*/

#ifndef _STREAM_BUFFER_HEADER
#define _STREAM_BUFFER_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>

typedef struct
{
  char *b;
  // begin is index of first char (unless begin >= end)
  // end is index of \0
  // size should be >= end+1 to allow for \0
  // (end-size) is the number of bytes in buffer
  size_t begin, end, size;
} CharBuffer;

// Buffer functions
#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) (0x1UL << (64 - __builtin_clzl(x)))
#endif

// Returns NULL if out of memory, @b otherwise
static inline CharBuffer* buffer_alloc(CharBuffer *b, size_t s)
{
  b->size = s < 16 ? 16 : s;
  if((b->b = malloc(b->size)) == NULL) return NULL;
  b->begin = b->end = 1;
  b->b[b->end] = b->b[b->size-1] = 0;
  return b;
}

static inline void buffer_dealloc(CharBuffer *b)
{
  free(b->b);
  memset(b, 0, sizeof(CharBuffer));
}

static inline CharBuffer* buffer_new(size_t s)
{
  CharBuffer *b = (CharBuffer*)malloc(sizeof(CharBuffer));
  if(b == NULL) return NULL;
  else if(buffer_alloc(b,s)) return b;
  else { free(b); return NULL; } /* couldn't malloc */
}

static inline void buffer_free(CharBuffer *cbuf)
{
  free(cbuf->b);
  free(cbuf);
}

// Resize a void pointer
// Adds one to len for NULL terminating byte
static inline void cbuffer_ensure_capacity(char **vbuf, size_t *sizeptr, size_t len)
{
  len++; // for nul byte
  if(*sizeptr < len) {
    *sizeptr = ROUNDUP2POW(len);
    if((*vbuf = realloc(*vbuf, *sizeptr)) == NULL) {
      fprintf(stderr, "Out of memory\n");
      exit(EXIT_FAILURE);
    }
  }
}

// len is the number of bytes you want to be able to store
// Adds one to len for NULL terminating byte
static inline void buffer_ensure_capacity(CharBuffer *cbuf, size_t len)
{
  cbuffer_ensure_capacity(&cbuf->b, &cbuf->size, len);
}

static inline void buffer_append_str(CharBuffer *buf, const char *str)
{
  size_t len = strlen(str);
  buffer_ensure_capacity(buf, buf->end+len);
  memcpy(buf->b+buf->end, str, len);
  buf->end += len;
  buf->b[buf->end] = 0;
}

static inline void buffer_append_char(CharBuffer *buf, char c)
{
  buffer_ensure_capacity(buf, buf->end+1);
  buf->b[buf->end++] = c;
  buf->b[buf->end] = '\0';
}

#define buffer_terminate(buf) ((buf)->b[(buf)->end] = 0)

#define buffer_reset(buf) do { \
  (buf)->begin = (buf)->end = 1; (buf)->b[1] = 0; \
} while(0)

#define buffer_len(buf) ((buf)->end - (buf)->begin)

static inline void buffer_chomp(CharBuffer *buf)
{
  while(buf->end > buf->begin &&
        (buf->b[buf->end-1] == '\n' || buf->b[buf->end-1] == '\r'))
  {
    buf->end--;
  }
  buf->b[buf->end] = 0;
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


// Define read for gzFile and FILE (unbuffered)
#define gzread2(gz,buf,len) gzread(gz,buf,(unsigned int)len)
#define fread2(f,buf,len) fread(buf,1,len,file)

#define gzgets2(gz,buf,len) gzgets(gz,buf,(int)(len))
#define fgets2(f,buf,len) fgets(buf,(int)(len),f)

// fgetc(f), gzgetc(gz) are already good to go
// fungetc(c,f), gzungetc(c,gz) are already good to go

#define fungetc(c,f) ungetc(c,f)

// Define readline for gzFile and FILE (unbuffered)
#define _func_readline(name,type_t,__gets) \
  static inline size_t name(type_t file, char **buf, size_t *len, size_t *size)\
  {                                                                            \
    if(*len+1 >= *size) cbuffer_ensure_capacity(buf, size, *len+1);            \
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
// both return -1 on error, 0 if nothing read, >0 on success
// offset of 1 so we can unget at least one char
// Beware: read-in buffer is not null-terminated
// Returns fail on error
// Can't have a buffer larger than 2^32 because that's all that gzFile can read
#define _READ_BUFFER(file,in,__read,fail) do                                   \
{                                                                              \
  size_t buf_size = (in)->size > UINT_MAX ? UINT_MAX : (in)->size;             \
  long _input = (long)__read(file,(in)->b+1,buf_size-1);                       \
  if(_input < 0) return fail;                                                  \
  (in)->end = 1+(size_t)_input; (in)->begin = 1;                               \
} while(0)

// Define getc for gzFile and FILE (buffered)
#define _func_getc_buf(fname,type_t,__read)                                    \
  static inline int fname(type_t file, CharBuffer *in)                         \
  {                                                                            \
    if(in->begin >= in->end) {                                                 \
      _READ_BUFFER(file,in,__read,-1);                                         \
      return in->end == in->begin ? -1 : in->b[in->begin++];                   \
    }                                                                          \
    return in->b[in->begin++];                                                 \
  }

_func_getc_buf(gzgetc_buf,gzFile,gzread2)
_func_getc_buf(fgetc_buf,FILE*,fread2)

// Define ungetc for buffers
// returns c if successful, otherwise -1
static inline int ungetc_buf(int c, CharBuffer *in)
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

#define _func_read_buf(fname,type_t,__read)                                    \
  static inline int fname(type_t file, void *ptr, size_t len, CharBuffer *in)  \
  {                                                                            \
    if(in->begin >= in->end) { _READ_BUFFER(file,in,__read,-1); }              \
    size_t remaining = len, next;                                              \
    while(in->end > in->begin && remaining > 0) {                              \
      next = in->end - in->begin;                                              \
      if(remaining <= next) next = remaining;                                  \
      memcpy(ptr, in->b+in->begin, next);                                      \
      in->begin += next; ptr = (char*)ptr + next;                              \
      if(remaining > next) { _READ_BUFFER(file,in,__read,-1); }                \
      remaining -= next;                                                       \
    }                                                                          \
    return (int)(len - remaining);                                             \
  }

_func_read_buf(gzread_buf,gzFile,gzread2)
_func_read_buf(fread_buf,FILE*,fread2)

// Define readline for gzFile and FILE (buffered)
#define _func_readline_buf(fname,type_t,__read)                                \
  static inline int fname(type_t file, CharBuffer *in,                         \
                          char **buf, size_t *len, size_t *size)               \
  {                                                                            \
    if(in->begin >= in->end) { _READ_BUFFER(file,in,__read,-1); }              \
    size_t offset, buffered, total_read = 0;                                   \
    while(in->end > in->begin)                                                 \
    {                                                                          \
      for(offset = in->begin; offset < in->end && in->b[offset++] != '\n'; ){} \
      buffered = offset - in->begin;                                           \
      cbuffer_ensure_capacity(buf, size, (*len)+buffered);                     \
      memcpy((*buf)+(*len), in->b+in->begin, buffered);                        \
      *len += buffered;                                                        \
      in->begin = offset;                                                      \
      total_read += buffered;                                                  \
      if((*buf)[*len-1] == '\n') break;                                        \
      _READ_BUFFER(file,in,__read,-1);                                         \
    }                                                                          \
    (*buf)[*len] = 0;                                                          \
    return (int)total_read;                                                    \
  }

_func_readline_buf(gzreadline_buf,gzFile,gzread2)
_func_readline_buf(freadline_buf,FILE*,fread2)

// Define buffered skipline
#define _func_skipline_buf(fname,ftype,__read)                                 \
  static inline int fname(ftype file, CharBuffer *in)                          \
  {                                                                            \
    if(in->begin >= in->end) { _READ_BUFFER(file,in,__read,-1); }              \
    size_t offset, skipped_bytes = 0;                                          \
    while(in->end > in->begin)                                                 \
    {                                                                          \
      for(offset = in->begin; offset < in->end && in->b[offset++] != '\n'; )   \
      skipped_bytes += offset - in->begin;                                     \
      in->begin = offset;                                                      \
      if(in->b[offset-1] == '\n') break;                                       \
      _READ_BUFFER(file,in,__read,-1);                                         \
    }                                                                          \
    return (int)skipped_bytes;                                                 \
  }

_func_skipline_buf(gzskipline_buf,gzFile,gzread2)
_func_skipline_buf(fskipline_buf,FILE*,fread2)

// Define buffered gzgets_buf, fgets_buf

// Reads upto len-1 bytes (or to the first \n if first) into str
// Adds null-terminating byte
#define _func_gets_buf(fname,type_t,__read) \
  static inline char* fname(type_t file, CharBuffer *in, char *str,            \
                            unsigned int len)                                  \
  {                                                                            \
    if(len == 0) return NULL;                                                  \
    if(len == 1) {str[0] = 0; return str; }                                    \
    if(in->begin >= in->end) { _READ_BUFFER(file,in,__read,NULL); }            \
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
      _READ_BUFFER(file,in,__read,NULL);                                       \
    }                                                                          \
    str[total_read] = 0;                                                       \
    return total_read == 0 ? NULL : str;                                       \
  }

_func_gets_buf(gzgets_buf,gzFile,gzread2)
_func_gets_buf(fgets_buf,FILE*,fread2)


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

#define fputs2(fh,c) fputs(c,fh)
#define gzputs2(gz,c) gzputs(gz,c)

#define fwrite2(fh,ptr,len) fwrite(ptr,len,1,fh)
#define gzwrite2(gz,ptr,len) gzwrite(gz,ptr,len)

/*
 Output (buffered)

// To do
fputc_buf(fh,buf,c)
gzputc_buf(gz,buf,c)
fputs_buf(fh,buf,str)
gzputs_buf(gz,buf,str)
fprintf_buf(fh,buf,fmt,...)
gzprintf_buf(gz,buf,fmt,...)
fwrite_buf(fh,buf,ptr,len)
gzwrite_buf(gz,buf,ptr,len)
buffer_flush(fh,buf)
buffer_gzflush(gz,buf)
*/


#endif
