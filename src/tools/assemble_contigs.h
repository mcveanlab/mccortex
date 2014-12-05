#ifndef ASSEMBLE_CONTIGS_H_
#define ASSEMBLE_CONTIGS_H_

#include "db_graph.h"
#include "contig_confidence.h"
#include "assemble_stats.h"

#include "seq_file/seq_file.h"

/**
 * Assemble contig for a given sample.
 *
 * @param seed_files If passed, use seed kmers from sequences. If not given,
 *                   iterate through the hash table.
 * @param contig_limit Stop after printing this many contigs, if zero no limit
 * @param visited Bit array to store visited nodes in. If not NULL, do not use a
 *                seed that has previously been visited. We do not clear this
 *                array.
 * @param seed_with_unused_paths If set, mark paths as used once entirely
 *                               contained in a contig. Unused paths are then
 *                               used to seed contigs.
 * @param min_confid Stop traversal if confidence drops below the given min
 *                   If less than 0 ignore.
 */
void assemble_contigs(size_t nthreads,
                      seq_file_t **seed_files, size_t num_seed_files,
                      size_t contig_limit, uint8_t *visited,
                      bool use_missing_info_check, bool seed_with_unused_paths,
                      double min_confid,
                      FILE *fout, const char *out_path,
                      AssembleContigStats *stats,
                      const ContigConfidenceTable *conf_table,
                      const dBGraph *db_graph, size_t colour);

#endif /* ASSEMBLE_CONTIGS_H_ */
