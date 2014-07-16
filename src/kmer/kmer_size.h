#ifndef KMER_SIZE_H_
#define KMER_SIZE_H_

// Files that are not compiled with MIN_KMER_SIZE and MAX_KMER_SIZE link to
// this object file and discover kmer size limits at run time

int get_min_kmer_size();
int get_max_kmer_size();

#endif /* KMER_SIZE_H_ */
