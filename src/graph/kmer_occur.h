#ifndef KMER_OCCUR_H_
#define KMER_OCCUR_H_

#include "seq_file.h"
#include "db_graph.h"

//
// This file provides a datastore for loading sequences and recording where
// each kmer occurs in the sequences. Used in breakpoint_caller.c.
//

// Limits on number of chromosomes and max chromosome length
// [30] max chroms    = 1,073,741,824
// [32] max chrom len = 4,294,967,296
#define KMER_OCCUR_MAX_CHROMS (1U<<30)
#define KMER_OCCUR_MAX_LEN (1UL<<32)

typedef struct {
  size_t id, length;
  const char *name;
} KOChrom;

// KOccur fits in one 64 bit word
typedef struct
{
  // Use a single bit to signify another entry for this kmer
  uint64_t next:1, orient:1, chrom:30, offset:32;
} KOccur;

// First we count the number of times a kmer occurs in the ref
// then we make a list
typedef union { uint64_t kcount; KOccur *first; } KONodeList;

typedef struct
{
  KOChrom *chroms; // Chromosomes in the reference genome
  KOccur *koccurs; // Kmer from ref that is in the graph
  KONodeList *klists; // one entry per hash entry
  size_t nchroms;
  char *chrom_name_buf;
} KOGraph;

typedef struct {
  uint64_t first, last; // 0-bases chromosome coordinates
  uint32_t qoffset, chrom; // qoffset some query offset
  // strand:
  //   0 => + (with ref, first <= last)
  //   1 => - (reverse ref, last <= first)
  bool strand, used; // used helps mark runs that are dropped
} KOccurRun;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(kmer_run_buf, KOccurRunBuffer, KOccurRun);

/**
 * Create a KOGraph from given sequence reads
 * BEWARE: We add the reads to the graph if add_missing_kmers is true
 * db_graph->col_edges can be NULL even if we are adding kmers
 * @param add_missing_kmers  If true, add kmers to the graph in colour ref_col
 **/
KOGraph kograph_create(const read_t *reads, size_t num_reads,
                       bool add_missing_kmers, size_t ref_col,
                       size_t num_threads, dBGraph *db_graph);

void kograph_free(KOGraph kograph);

// Get KOccur* to first occurance of a kmer in sequence
#define kograph_get(kograph,hkey) ((kograph).klists[hkey].first)

#define kograph_occurs(kograph,hkey) ((kograph).klists[hkey].first != NULL)

// Get the chromosome from which a kmer came (occur can be KOccurRun or KOccur)
#define kograph_chrom(kograph,occur) ((kograph).chroms[(occur).chrom])

#define korun_len(run) ((run).strand == STRAND_PLUS ? (run).last-(run).first+1 \
                                                    : (run).first-(run).last+1)

void koruns_sort_by_qoffset(KOccurRun *runs, size_t n);

// Filter regions down to only those that stretch the whole distance
// Does not reset either korun or runs_ended - only adds to runs_ended
// korun can only get shorter as KOccurRuns finish
// `qoffset` is used for offset of new runs starting at nodes[0]
void kograph_filter_extend(KOGraph kograph,
                           const dBNode *nodes, size_t num_nodes, bool forward,
                           size_t min_len, size_t qoffset,
                           KOccurRunBuffer *korun,
                           KOccurRunBuffer *runs_ended,
                           bool pickup_at_first_node);

// Mostly used for debugging
void korun_print(KOccurRun run, size_t kmer_size, FILE *fout);

// Get string representation of multiple runs, comma separated
//   e.g. "chromid:1:17-5:-, chromid:1:37-47:+"
// Does not print new line
// Mostly used for debugging
void koruns_print(KOccurRun *run, size_t n, size_t kmer_size, FILE *fout);

void korun_gzprint(gzFile gzout, size_t kmer_size,
                   KOGraph kograph, KOccurRun korun,
                   size_t first_kmer_idx, size_t kmer_offset);

void koruns_gzprint(gzFile gzout, size_t kmer_size, KOGraph kograph,
                    const KOccurRun *koruns, size_t n,
                    size_t first_kmer_idx, size_t kmer_offset);

// src, dst can point to the same place
// returns number of elements added
static inline
size_t koruns_filter(KOccurRun *src, size_t n, KOccurRun *dst, size_t min_kmers)
{
  size_t i, j;
  for(i = j = 0; i < n; i++)
    if(korun_len(src[i]) >= min_kmers)
      dst[j++] = src[i];

  return j;
}

static inline void koruns_reverse(KOccurRun *src, size_t n, size_t nkmers)
{
  size_t i, j, len;
  for(i = j = 0; i < n; i++) {
    len = (src[i].strand == STRAND_PLUS ? src[i].last - src[i].first
                                        : src[i].first - src[i].last);
    SWAP(src[i].first, src[i].last);
    src[i].strand = !src[i].strand;
    src[i].qoffset = nkmers-1-src[i].qoffset-len;
  }
}

#endif /* KMER_OCCUR_H_ */
