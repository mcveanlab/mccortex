#ifndef GPATH_SAVE_H_
#define GPATH_SAVE_H_

#include "db_graph.h"
#include "gpath_subset.h"
#include "cJSON/cJSON.h"

/*
// File format:
<JSON_HEADER>
kmer [num] .. ignored
[FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. ignored
*/

// @hdrs is array of JSON headers of input files
void gpath_save(gzFile gzout, const char *path, size_t nthreads,
                cJSON **hdrs, size_t nhdrs,
                const ZeroSizeBuffer *contig_hists, size_t ncols,
                dBGraph *db_graph);

// Save paths for a single kmer
// @subset is temporary memory
void gpath_fwrite_single_kmer(hkey_t hkey, FILE *fout,
                              GPathSubset *subset,
                              const dBGraph *db_graph);

#endif /* GPATH_SAVE_H_ */
