#include "global.h"
#include "seqout.h"
#include "file_util.h"

#include <libgen.h> // basename(), dirname()

// Malloc and return path with given suffix
static char* _seqout_alloc_path(char *out_base, const char *suffix)
{
  size_t len1 = strlen(out_base), len2 = strlen(suffix);
  char *path = ctx_malloc((len1+len2+1) * sizeof(char));
  memcpy(path, out_base, len1);
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
bool seqout_open(SeqOutput *seqout, char *out_base, bool use_fq, bool is_pe)
{
  memset(seqout, 0, sizeof(SeqOutput));

  seqout->use_fq = use_fq;
  seqout->is_pe = is_pe;
  seqout->path_se = _seqout_alloc_path(out_base, use_fq ? ".fq.gz" : ".fa.gz");
  if((seqout->gzout_se = _seqout_open(seqout->path_se)) == NULL) return false;

  if(is_pe) {
    seqout->path_pe[0] = _seqout_alloc_path(out_base, use_fq ? ".1.fq.gz" : ".1.fa.gz");
    seqout->path_pe[1] = _seqout_alloc_path(out_base, use_fq ? ".2.fq.gz" : ".2.fa.gz");
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

#define print_read(r,use_fq,gzout)                \
do {                                              \
  if(use_fq) seq_gzprint_fastq(r, gzout, 0);      \
  else       seq_gzprint_fasta(r, gzout, 0);      \
} while(0)

void seqout_print(SeqOutput *seqout, const read_t *r1, const read_t *r2)
{
  if(r2 == NULL) {
    pthread_mutex_lock(&seqout->lock_se);
    print_read(r1, seqout->use_fq, seqout->gzout_se);
    pthread_mutex_unlock(&seqout->lock_se);
  } else {
    pthread_mutex_lock(&seqout->lock_pe);
    print_read(r1, seqout->use_fq, seqout->gzout_pe[0]);
    print_read(r2, seqout->use_fq, seqout->gzout_pe[1]);
    pthread_mutex_unlock(&seqout->lock_pe);
  }
}
