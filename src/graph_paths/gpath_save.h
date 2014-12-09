#ifndef GPATH_SAVE_H_
#define GPATH_SAVE_H_

#include "db_graph.h"
#include "db_node.h"
#include "gpath_subset.h"
#include "cJSON/cJSON.h"

/*
// File format:
<JSON_HEADER>
kmer [num] .. ignored
[FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. ignored
*/

extern const char ctp_explanation_comment[];

cJSON* gpath_save_mkhdr(const char *path,
                        cJSON **hdrs, size_t nhdrs,
                        const ZeroSizeBuffer *contig_hists, size_t ncols,
                        const dBGraph *db_graph);

/**
 * Print paths to a string buffer. Paths are sorted before being written.
 *
 * @param hkey    All paths associated with hkey are written to the buffer
 * @param sbuf    paths are written this string buffer
 * @param subset  is a temp variable that is reused each time
 * @param nbuf    temporary buffer, if not NULL, used to add seq=... to output
 * @param jposbuf temporary buffer, if not NULL, used to add juncpos=... to output
 */
void gpath_save_sbuf(hkey_t hkey, StrBuf *sbuf, GPathSubset *subset,
                     dBNodeBuffer *nbuf, SizeBuffer *jposbuf,
                     const dBGraph *db_graph);

/**
 * Save paths to a file.
 * @param hdrs is array of JSON headers of input files
 */
void gpath_save(gzFile gzout, const char *path,
                size_t nthreads, bool save_path_seq,
                cJSON **hdrs, size_t nhdrs,
                const ZeroSizeBuffer *contig_hists, size_t ncols,
                dBGraph *db_graph);

#endif /* GPATH_SAVE_H_ */
