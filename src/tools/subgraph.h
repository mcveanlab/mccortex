#ifndef SUBGRAPH_H_
#define SUBGRAPH_H_

#include "seq_file.h"
#include "db_graph.h"

//
// Breadth first search from seed kmers,
// then remove all kmers that weren't touched
//

// `dist` is how many steps away from seed kmers to take
// `invert`, if true, means only save kmers not touched
// `fringe_mem` is how many bytes can be used to remember the fringe of
// the breadth first search (we use 8 bytes per kmer)
// `kmer_mask` should be a bit array (one bit per kmer) of zero'd memory
void subgraph_from_reads(dBGraph *db_graph, size_t dist,
                         bool invert, bool grab_supernodes,
                         size_t fringe_mem, uint64_t *kmer_mask,
                         seq_file_t **files, size_t num_files);

// `dist` is how many steps away from seed kmers to take
// `invert`, if true, means only save kmers not touched
// `fringe_mem` is how many bytes can be used to remember the fringe of
// the breadth first search (we use 8 bytes per kmer)
// `kmer_mask` should be a bit array (one bit per kmer) of zero'd memory
void subgraph_from_seq(dBGraph *db_graph, size_t dist,
                       bool invert, bool grab_supernodes,
                       size_t fringe_mem, uint64_t *kmer_mask,
                       char **seqs, size_t *seqlens, size_t num_seqs);

#endif /* SUBGRAPH_H_ */
