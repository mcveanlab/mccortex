#ifndef SEQ_OUTPUT_H_
#define SEQ_OUTPUT_H_

//
// FASTA/Q output
//
// writes to files:
//   <out>.fa.gz
//   <out>.1.fa.gz
//   <out>.2.fa.gz
//

#include "seq_file.h"

typedef struct {
  char *path_se, *path_pe[2];
  gzFile gzout_se, gzout_pe[2];
  pthread_mutex_t lock_se, lock_pe;
  bool is_pe; // if we have X.{1,2}.fq.gz as well as X.fq.gz
  seq_format fmt; // output format
} SeqOutput;

// Returns true on success, false on failure
// fmt may be: SEQ_FMT_FASTQ, SEQ_FMT_FASTA, SEQ_FMT_PLAIN
bool seqout_open(SeqOutput *seqout, char *out_base, seq_format fmt, bool is_pe);

// Free memory
// @rm if true, delete files as well
void seqout_close(SeqOutput *output, bool rm);

void seqout_print(SeqOutput *output, const read_t *r1, const read_t *r2);

#endif /* SEQ_OUTPUT_H_ */
