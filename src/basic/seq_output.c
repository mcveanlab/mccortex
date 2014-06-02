#include "global.h"
#include "seq_output.h"
#include "file_util.h"

//
// FASTA output
//

void seq_output_alloc(SeqOutput *out)
{
  strbuf_alloc(&out->path_se, 512);
  strbuf_alloc(&out->path_pe[0], 512);
  strbuf_alloc(&out->path_pe[1], 512);
  out->gzout_se = out->gzout_pe[0] = out->gzout_pe[1] = NULL;
  out->output_pe = false;
  if(pthread_mutex_init(&out->lock_se, NULL) != 0) die("Mutex init failed");
  if(pthread_mutex_init(&out->lock_pe, NULL) != 0) die("Mutex init failed");
}

void seq_output_dealloc(SeqOutput *out)
{
  if(out->gzout_se) gzclose(out->gzout_se);
  if(out->gzout_pe[0]) gzclose(out->gzout_pe[0]);
  if(out->gzout_pe[1]) gzclose(out->gzout_pe[1]);
  out->gzout_se = out->gzout_pe[0] = out->gzout_pe[1] = NULL;
  strbuf_dealloc(&out->path_se);
  strbuf_dealloc(&out->path_pe[0]);
  strbuf_dealloc(&out->path_pe[1]);
  out->output_pe = false;
  pthread_mutex_destroy(&out->lock_se);
  pthread_mutex_destroy(&out->lock_pe);
}

void seq_output_set_paths(SeqOutput *out, const char *base, bool pe)
{
  strbuf_reset(&out->path_se);
  strbuf_reset(&out->path_pe[0]);
  strbuf_reset(&out->path_pe[1]);

  // Single ended
  strbuf_set(&out->path_se, base);
  strbuf_append_str(&out->path_se, ".fa.gz");

  if(pe) {
    strbuf_set(&out->path_pe[0], base);
    strbuf_set(&out->path_pe[1], base);
    strbuf_append_str(&out->path_pe[0], ".1.fa.gz");
    strbuf_append_str(&out->path_pe[1], ".2.fa.gz");
  }

  out->output_pe = pe;
}

// Close and delete opened files
void seq_output_delete(SeqOutput *out)
{
  #define test_and_delete(gz,pathbuf) do {                     \
    if((gz) != NULL) { gzclose(gz); unlink((pathbuf)->buff); } \
  } while(0)

  test_and_delete(out->gzout_se, &out->path_se);
  test_and_delete(out->gzout_pe[0], &out->path_pe[0]);
  test_and_delete(out->gzout_pe[1], &out->path_pe[1]);
  out->gzout_se = out->gzout_pe[0] = out->gzout_pe[1] = NULL;

  #undef test_and_delete
}

// Call warn() with an error if any of the proposed output files already exist
// Returns:
//  - true if one or more files already exist
//  - false if no files already exist
bool seq_output_files_exist_check(const SeqOutput *out)
{
  size_t i;
  bool exists = false;

  if(futil_file_exists(out->path_se.buff)) {
    warn("SE file already exists: %s", out->path_se.buff);
    exists = true;
  }

  if(out->output_pe) {
    for(i = 0; i < 2; i++) {
      if(futil_file_exists(out->path_pe[i].buff)) {
        warn("PE file already exists: %s", out->path_pe[i].buff);
        exists = true;
      }
    }
  }

  return exists;
}

// Returns true on success, false on error
bool seq_output_open(SeqOutput *out)
{
  size_t i;

  out->gzout_se = gzopen(out->path_se.buff, "w");
  if(out->gzout_se == NULL) {
    warn("Cannot write to: %s", out->path_se.buff);
    return false;
  }

  if(out->output_pe) {
    for(i = 0; i < 2; i++) {
      out->gzout_pe[i] = gzopen(out->path_pe[i].buff, "w");
      if(out->gzout_pe[i] == NULL) {
        warn("Cannot write to: %s", out->path_pe[i].buff);
        seq_output_delete(out); // remove existing files
        return false;
      }
    }
  }

  return true;
}
