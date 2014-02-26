#ifndef DB_ALIGNMENT_H_
#define DB_ALIGNMENT_H_

#include "db_graph.h"
#include "db_node.h"
#include "seq_file.h"

// Buffer of uint32_t
#include "objbuf_macro.h"
create_objbuf(uint32_buf,uint32Buffer,uint32_t);

typedef struct
{
  dBNodeBuffer nodes; // Nodes found in the graph
  uint32Buffer gaps; // gap (num nodes) preceeding each node covered by sequence
  size_t r1enderr, r2enderr; // how many nodes are lost at the end of reads 1,2
  // index of buf where r2 starts
  // (-1 if no r2, nodes.len if r2 did not provide any nodes)
  // int r2strtidx;
  // 0 if no r1 nodes, nodes.len if
  size_t r2strtidx;
  // whether we were passed a second read
  boolean passed_r2;
  // used_r1, used_r2 are whether or not we got nodes from each read
  boolean used_r1, used_r2;
  // whether or not there are sequencing gaps
  boolean seq_gaps;
  // if used_r1 && used_r2
  // gap between r1 last nodes[r2strtindx-1] .. nodes[r2strtindx]
  // = r1enderr + insgapsize + gaps[r2strtindx]
} dBAlignment;

// Estimate memory required
size_t db_alignment_est_mem();

void db_alignment_alloc(dBAlignment *alignment);
void db_alignment_dealloc(dBAlignment *alignment);

void db_alignment_from_reads(dBAlignment *alignment,
                             const read_t *r1, const read_t *r2,
                             uint8_t qcutoff1, uint8_t qcutoff2,
                             uint8_t hp_cutoff, const dBGraph *db_graph);

// Returns index of node just after next gap,
// or aln->nodes.len if no more gaps
size_t db_alignment_next_gap(const dBAlignment *aln, size_t start);

//
// Debugging
//

void db_alignment_print(const dBAlignment *aln, const dBGraph *db_graph);

// Check all edges between ungapped adjacent nodes
boolean db_alignment_check_edges(const dBAlignment *aln, const dBGraph *graph);

// dBKmer stores redundant data from the graph to speed up processing
// typedef struct
// {
//   dBNode node;
//   BinaryKmer bkmer;
//   Edges edges;
// } dBKmer;

// #include "objbuf_macro.h"
// create_objbuf(db_kmer_buf,dBKmerBuffer,dBKmer)

#endif /* DB_ALIGNMENT_H_ */
