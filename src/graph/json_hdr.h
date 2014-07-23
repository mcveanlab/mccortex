#ifndef JSON_HDR_H_
#define JSON_HDR_H_

#include "db_graph.h"
#include "cJSON/cJSON.h"

#define MAX_JSON_HDR_BYTES (1<<20) /* 1M max json header */

// Read json header into hdrstr
// read from FILE* or gzFile, whichever is not NULL
void json_hdr_read(FILE *fh, gzFile gz, const char *path, StrBuf *hdrstr);

// Add standard header fields to a json header
// Merge commands from input files @hdrs
// @param path is the path of the file we are writing to
void json_hdr_add_std(cJSON *json, const char *path,
                      cJSON **hdrs, size_t nhdrs,
                      const dBGraph *db_graph);

void json_hdr_gzprint(cJSON *json, gzFile gzout);

#endif /* JSON_HDR_H_ */
