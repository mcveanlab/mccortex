#ifndef KMER_OCCUR_H_
#define KMER_OCCUR_H_

#include "seq_file.h"
#include "db_graph.h"

//
// This file provides a datastore for loading sequences and recording where
// each kmer occurs in the sequences. Used in breakpoint_caller.c.
//

typedef struct {
  size_t id, length;
  const char *name;
} KOChrom;

// 8+4=12 bytes
struct KONodeListStruct {
  uint64_t start;
  uint32_t count;
} __attribute__((packed));

typedef struct KONodeListStruct KONodeList;

typedef struct
{
  uint64_t orient:1, chrom:23, offset:40;
} KOccur;

typedef struct
{
  KOChrom *chroms; // Chromosomes in the reference genome
  KOccur *koccurs; // Kmer from ref that is in the graph
  KONodeList *klists; // one entry per hash entry
  size_t nchroms;
  char *chrom_name_buf;
} KOGraph;

#include "objbuf_macro.h"
create_objbuf(kmer_occur_buf, KOccurBuffer, KOccur);

// Get KOccur* to first occurance of a kmer in sequence
#define kograph_get(kograph,hkey) ((kograph).koccurs + (kograph).klists[hkey].start)

// Get the number of times a kmer is seen in the sequence
#define kograph_num(kograph,hkey) ((kograph).klists[hkey].count)

#define kograph_get_check(kograph,hkey) \
        (!kograph_num(kograph,hkey) ? NULL : kograph_get(kograph,hkey))

// Get the chromosome from which a kmer came
#define kograph_chrom(kograph,occur) (&(kograph).chroms[(occur)->chrom])

// We add the reads to the graph if `add_missing_kmers` is true
KOGraph kograph_create(const read_t *reads, size_t num_reads,
                       bool add_missing_kmers, size_t num_threads,
                       dBGraph *db_graph);

void kograph_free(KOGraph kograph);

// occur0 should be before occur1
// kolist0 and kolist1 should already be sorted
void kograph_get_kmer_run(const KOccur *restrict kolist0, size_t num0,
                          const KOccur *restrict kolist1, size_t num1,
                          size_t dist, KOccurBuffer *restrict kobuf);

#endif /* KMER_OCCUR_H_ */
