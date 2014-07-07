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
  bool use_fq; // .fa.gz or .fq.gz
} SeqOutput;

// Returns true on success, false on failure
bool seqout_open(SeqOutput *seqout, char *out_base, bool use_fq, bool is_pe);

// Free memory
// @rm if true, delete files as well
void seqout_close(SeqOutput *output, bool rm);

void seqout_print(SeqOutput *output, const read_t *r1, const read_t *r2);

#endif /* SEQ_OUTPUT_H_ */
