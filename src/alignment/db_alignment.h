#ifndef DB_ALIGNMENT_H_
#define DB_ALIGNMENT_H_

#include "db_graph.h"
#include "db_node.h"
#include "common_buffers.h" // Buffer of uint32_t

#include "seq_file/seq_file.h"

typedef struct
{
  size_t r1bases, r2bases; // number of bases in input reads
  dBNodeBuffer nodes; // Nodes found in the graph
  Int32Buffer rpos; // Position of kmer start in read e.g {0,1,5,6,0,1,2} in PE
  size_t r1enderr, r2enderr; // how many nodes are lost at the end of reads 1,2
  // index of buf where r2 starts
  // 0 if no r1 nodes, nodes.len if no r2 does not provide any nodes
  size_t r2strtidx;
  // whether we were passed a second read
  bool passed_r2;
  // used_r1, used_r2 are whether or not we got nodes from each read
  bool used_r1, used_r2;
  // whether or not there are sequencing gaps
  bool seq_gaps;
  // if used_r1 && used_r2
  // gap between r1 and r2: nodes[r2strtidx-1] .. nodes[r2strtidx]
  // = r1enderr + insgapsize + rpos[r2strtidx]
  int colour; // -1 if colour agnostic, otherwise only nodes in colour used
} dBAlignment;

// Estimate memory required
size_t db_alignment_est_mem();

void db_alignment_alloc(dBAlignment *alignment);
void db_alignment_dealloc(dBAlignment *alignment);

#define db_aln_r1enderr(aln,k) ((aln)->r1bases - \
  ((aln)->r2strtidx > 0 ? (aln)->rpos.b[(aln)->r2strtidx-1] + (k) : 0))

// if colour is -1 aligns to all colours, otherwise aligns to given colour only
// Assumes both reads are in FF orientation
void db_alignment_from_reads(dBAlignment *alignment,
                             const read_t *r1, const read_t *r2,
                             uint8_t qcutoff1, uint8_t qcutoff2,
                             uint8_t hp_cutoff,
                             const dBGraph *db_graph, int colour);

/*
 * Get position after current segment. A segement is a stretch of kmers aligned
 * to the graph with no gaps.
 * @param aln           alignment of read to the graph
 * @param start         offset in the alignment that starts this segments
 * @param missing_edge  set to 1 if the segment was ended due to a missing edge
 * @return index of node after next gap or aln->nodes.len if no more gaps
 */
size_t db_alignment_next_gap(const dBAlignment *aln, size_t start,
                             bool *missing_edge, const dBGraph *db_graph);

// @return true iff alignment has no gaps, all kmers found
bool db_alignment_is_perfect(const dBAlignment *aln);

//
// Debugging
//

void db_alignment_print(const dBAlignment *aln);

// Check all edges between ungapped adjacent nodes
// bool db_alignment_check_edges(const dBAlignment *aln, const dBGraph *graph);

// dBKmer stores redundant data from the graph to speed up processing
// typedef struct
// {
//   dBNode node;
//   BinaryKmer bkmer;
//   Edges edges;
// } dBKmer;

// #include "madcrowlib/madcrow_buffer.h"
// madcrow_buffer(db_kmer_buf,dBKmerBuffer,dBKmer)

#endif /* DB_ALIGNMENT_H_ */
