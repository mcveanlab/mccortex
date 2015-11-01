#ifndef BREAKPOINT_CALLER_H_
#define BREAKPOINT_CALLER_H_

#include "db_graph.h"

#include "seq_file/seq_file.h"
#include "cJSON/cJSON.h"

#define BREAKPOINT_FORMAT_VERSION 3

// Require 5 kmers on the reference before and after breakpoint
#define DEFAULT_MIN_REF_NKMERS 5
// By default stop following the ref after 1000 kmers
#define DEFAULT_MAX_REF_NKMERS 1000

/**
 * Load reference kmers into the graph, make breakpoint calls and write out.
 *
 * @param nthreads      number of threads to use
 * @param ref_col       colour to load reference sequence into
 * @param gzout         gzFile to print breakpoints to
 * @param out_path      path to output file that gzout points to
 * @param reads         reference sequence to load into the graph
 * @param num_reads     number of reference contigs
 * @param seq_paths     paths to the files which the ref reads where loaded from
 * @param num_seq_paths number of seq_paths
 * @param load_ref_edges whether or not to load edges from the ref
 * @param min_ref_flank num of kmers required to flank breakpoint on ref
 * @param hdrs          JSON headers of input files
 * @param nhdrs         number of JSON headers in hdrs
 * @param db_graph      de Bruijn graph to use
 **/
void breakpoints_call(size_t nthreads, size_t ref_col,
                      gzFile gzout, const char *out_path,
                      const read_t *reads, size_t num_reads,
                      char **seq_paths, size_t num_seq_paths,
                      bool load_ref_edges,
                      size_t min_ref_flank, size_t max_ref_flank,
                      cJSON **hdrs, size_t nhdrs,
                      dBGraph *db_graph);

#endif /* BREAKPOINT_CALLER_H_ */
