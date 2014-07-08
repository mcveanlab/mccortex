#include "global.h"
#include "seqout.h"
#include "file_util.h"

#include <libgen.h> // basename(), dirname()

// Malloc and return path with given suffix
// @pe 0 if se, 1/2 if one of a pair (out.fq.gz out.1.fq.gz, out.2.fq.gz)
static char* _seqout_alloc_path(char *out_base, int pe, const char *suffix)
{
  ctx_assert2(pe >= 0 && pe <= 2, "bad PE param: %i", pe);
  size_t len1 = strlen(out_base), len2 = strlen(suffix);
  char *path = ctx_malloc((len1+2+len2+1) * sizeof(char)); // +2 for PE '.1'
  memcpy(path, out_base, len1);
  if(pe) { path[len1++] = '.'; path[len1++] = '0'+pe; }
  memcpy(path+len1, suffix, len2);
  path[len1+len2] = '\0';
  return path;
}

// Returns path to file or NULL if file already exists and futil_get_force()
// returns false
// Creates directories as required
static gzFile _seqout_open(const char *path)
{
  gzFile gzout;

  if(!futil_get_force() && futil_file_exists(path)) {
    warn("Output file already exists: %s", path);
    return NULL;
  }

  // dirname, basename may modify string, so make copy
  char *pathcpy = strdup(path);
  char *fname = basename(pathcpy);

  if(path[0] == '\0' || path[strlen(path)-1] == '\0' ||
     fname[0] == '/' || fname[0] == '.')
  {
    warn("Bad seqout name: %s", path);
    free(pathcpy);
    return NULL;
  }

  strcpy(pathcpy, path);
  char *dir = dirname(pathcpy);
  futil_mkpath(dir, 0777);
  free(pathcpy);

  if((gzout = gzopen(path, "w")) == NULL) {
    warn("Cannot open %s", path);
    return NULL;
  }

  // Set buffer size
  #if ZLIB_VERNUM >= 0x1240
    gzbuffer(gzout, DEFAULT_IO_BUFSIZE);
  #endif

  return gzout;
}

// Returns true on success, false on failure
bool seqout_open(SeqOutput *seqout, char *out_base, seq_format fmt, bool is_pe)
{
  memset(seqout, 0, sizeof(SeqOutput));

  seqout->fmt = fmt;
  seqout->is_pe = is_pe;
  const char *ext = NULL;

  switch(fmt) {
    case SEQ_FMT_FASTQ: ext = ".fq.gz";  break;
    case SEQ_FMT_FASTA: ext = ".fa.gz";  break;
    case SEQ_FMT_PLAIN: ext = ".txt.gz"; break;
    default: die("Invalid format: %i", (int)fmt);
  }

  seqout->path_se = _seqout_alloc_path(out_base, 0, ext);
  if((seqout->gzout_se = _seqout_open(seqout->path_se)) == NULL) return false;

  if(is_pe) {
    seqout->path_pe[0] = _seqout_alloc_path(out_base, 1, ext);
    seqout->path_pe[1] = _seqout_alloc_path(out_base, 2, ext);
    if((seqout->gzout_pe[0] = _seqout_open(seqout->path_pe[0])) == NULL) return false;
    if((seqout->gzout_pe[1] = _seqout_open(seqout->path_pe[1])) == NULL) return false;
  }

  if(pthread_mutex_init(&seqout->lock_se, NULL) != 0) die("Mutex init failed");
  if(pthread_mutex_init(&seqout->lock_pe, NULL) != 0) die("Mutex init failed");

  return true;
}

// Free memory
// @rm if true, delete files as well
void seqout_close(SeqOutput *seqout, bool rm)
{
  // Clean up seqout
  if(seqout->gzout_se != NULL) { gzclose(seqout->gzout_se); }
  if(seqout->gzout_pe[0] != NULL) { gzclose(seqout->gzout_pe[0]); }
  if(seqout->gzout_pe[1] != NULL) { gzclose(seqout->gzout_pe[1]); }
  if(rm) {
    if(seqout->gzout_se != NULL && unlink(seqout->path_se) != 0)
      warn("Cannot delete file %s", seqout->path_se);
    if(seqout->gzout_pe[0] != NULL && unlink(seqout->path_pe[0]) != 0)
      warn("Cannot delete file %s", seqout->path_pe[0]);
    if(seqout->gzout_pe[1] != NULL && unlink(seqout->path_pe[1]) != 0)
      warn("Cannot delete file %s", seqout->path_pe[1]);
  }
  ctx_free(seqout->path_se);
  ctx_free(seqout->path_pe[0]);
  ctx_free(seqout->path_pe[1]);
  pthread_mutex_destroy(&seqout->lock_se);
  pthread_mutex_destroy(&seqout->lock_pe);
  memset(seqout, 0, sizeof(SeqOutput));
}

#define print_read(r,fmt,gzout)                                                \
do {                                                                           \
  switch(fmt) {                                                                \
    case SEQ_FMT_PLAIN: gzputs(gzout, r->seq.b); gzputc(gzout, '\n'); break;   \
    case SEQ_FMT_FASTA: seq_gzprint_fasta(r, gzout, 0); break;                 \
    case SEQ_FMT_FASTQ: seq_gzprint_fastq(r, gzout, 0); break;                 \
    default: die("Invalid output format: %i", fmt);                            \
  }                                                                            \
} while(0)

void seqout_print(SeqOutput *seqout, const read_t *r1, const read_t *r2)
{
  if(r2 == NULL) {
    pthread_mutex_lock(&seqout->lock_se);
    print_read(r1, seqout->fmt, seqout->gzout_se);
    pthread_mutex_unlock(&seqout->lock_se);
  } else {
    pthread_mutex_lock(&seqout->lock_pe);
    print_read(r1, seqout->fmt, seqout->gzout_pe[0]);
    print_read(r2, seqout->fmt, seqout->gzout_pe[1]);
    pthread_mutex_unlock(&seqout->lock_pe);
  }
}
