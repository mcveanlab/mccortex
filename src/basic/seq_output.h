#ifndef SEQ_OUTPUT_H_
#define SEQ_OUTPUT_H_

//
// FASTA output
//
// writes to files:
//   <out>.fa.gz
//   <out>.1.fa.gz
//   <out>.2.fa.gz
//

typedef struct {
  StrBuf path_se, path_pe[2];
  gzFile gzout_se, gzout_pe[2];
  pthread_mutex_t lock_se, lock_pe;
  bool output_pe; // false => only one output file; true => three output files
} SeqOutput;

void seq_output_alloc(SeqOutput *out);

// Close files, free memory
void seq_output_dealloc(SeqOutput *out);

// Set paths to <base>.fa.gz <base>.1.fa.gz <base>.2.fa.gz
void seq_output_set_paths(SeqOutput *out, const char *base, bool pe);

// Close and delete opened files
void seq_output_delete(SeqOutput *out);

// Call warn() with an error if any of the proposed output files already exist
// Returns:
//  - true if one or more files already exist
//  - false if no files already exist
bool seq_output_files_exist_check(const SeqOutput *out);

// Returns true on success, false on error
// If cannot open file, removes opened files
bool seq_output_open(SeqOutput *out);

#endif /* SEQ_OUTPUT_H_ */
