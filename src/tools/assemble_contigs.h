#ifndef ASSEMBLE_CONTIGS_H_
#define ASSEMBLE_CONTIGS_H_

#include "db_graph.h"
#include "contig_confidence.h"
#include "assemble_stats.h"

#include "seq_file/seq_file.h"

/*!
  @param contig_limit Stop after printing this many contigs, if zero no limit
  @param visited Bit array to store visited nodes in. If not NULL, do not use a
                 seed that has previously been visited. We do not clear this
                 array.
 */
void assemble_contigs(size_t nthreads,
                      seq_file_t **seed_files, size_t num_seed_files,
                      size_t contig_limit, uint8_t *visited,
                      FILE *fout, const char *out_path,
                      AssembleContigStats *stats,
                      const ContigConfidenceTable *conf_table,
                      const dBGraph *db_graph, size_t colour);

#endif /* ASSEMBLE_CONTIGS_H_ */
