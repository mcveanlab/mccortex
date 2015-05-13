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
[FR] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. ignored
*/

extern const char ctp_explanation_comment[];

/**
 * Generate a JSON header object for a .ctp file
 * @param path        path to output file
 * @param cmdstr      name of the command being run, to be used to add @cmdhdr
 * @param cmdhdr      JSON header to add under current command->@cmdstr
 *                    If cmdstr and cmdhdr are both NULL they are ignored
 * @param hdrs        array of JSON headers of input files
 * @param nhdrs       number of elements in @hdrs
 * @param contig_hist histgram of read contig lengths
 * @param hist_len    length of array contig_hist
 */
cJSON* gpath_save_mkhdr(const char *path,
                        const char *cmdstr, cJSON *cmdhdr,
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
 * @param cmdstr  name of the command being run, to be used to add @cmdhdr
 * @param cmdhdr  JSON header to add under current command->@cmdstr
 *                If cmdstr and cmdhdr are both NULL they are ignored
 * @param hdrs    array of JSON headers of input files
 * @param nhdrs   number of elements in @hdrs
 */
void gpath_save(gzFile gzout, const char *path,
                size_t nthreads, bool save_path_seq,
                const char *cmdstr, cJSON *cmdhdr,
                cJSON **hdrs, size_t nhdrs,
                const ZeroSizeBuffer *contig_hists, size_t ncols,
                dBGraph *db_graph);

#endif /* GPATH_SAVE_H_ */
