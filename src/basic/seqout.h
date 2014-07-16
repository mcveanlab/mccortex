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
// file extensions are: <O>.fq.gz, <O>.fa.gz, <O>.txt.gz
bool seqout_open(SeqOutput *seqout, char *out_base, seq_format fmt, bool is_pe);

// Free memory
// @rm if true, delete files as well
void seqout_close(SeqOutput *output, bool rm);

void seqout_print(SeqOutput *output, const read_t *r1, const read_t *r2);

static inline int seqout_str2fmt(const char *str)
{
  ctx_assert(str != NULL);
  if(!strcasecmp(str,"fq") || !strcasecmp(str,"fastq")) return SEQ_FMT_FASTQ;
  if(!strcasecmp(str,"fa") || !strcasecmp(str,"fasta")) return SEQ_FMT_FASTA;
  if(!strcasecmp(str,"plain") || !strcasecmp(str,"txt")) return SEQ_FMT_PLAIN;
  if(!strcasecmp(str,"sam")) return SEQ_FMT_SAM;
  if(!strcasecmp(str,"bam")) return SEQ_FMT_BAM;
  return -1;
}

// str must be at least 6 bytes long
static inline char* seqout_fmt2str(seq_format fmt, char *str)
{
  ctx_assert(str != NULL);
  switch(fmt) {
    case SEQ_FMT_FASTQ: strcpy(str, "FASTQ"); break;
    case SEQ_FMT_FASTA: strcpy(str, "FASTA"); break;
    case SEQ_FMT_PLAIN: strcpy(str, "PLAIN"); break;
    case SEQ_FMT_SAM:   strcpy(str, "SAM"); break;
    case SEQ_FMT_BAM:   strcpy(str, "BAM"); break;
    default: die("Bad format: %i", fmt);
  }
  return str;
}

#endif /* SEQ_OUTPUT_H_ */
