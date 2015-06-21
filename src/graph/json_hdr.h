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
void json_hdr_make_std(cJSON *json, const char *path,
                       cJSON **hdrs, size_t nhdrs,
                       const dBGraph *db_graph);

// Add current command to a header
void json_hdr_add_curr_cmd(cJSON *json, const char *path);

/**
 * Add tags to current command under @field
 * @param json JSON header to add to. Must already have current command added -
 *             call json_hdr_add_curr_cmd() to add it first.
 * @param field name of field
 * @param add JSON objects to add to JSON object field
 * @param nadd number of @add objects
 */
void json_hdr_augment_cmd(cJSON *json, const char *cmdstr,
                          const char *field, cJSON *add);

void json_hdr_gzprint(cJSON *json, gzFile gzout);
void json_hdr_fprint(cJSON *json, FILE *fout);

// Get values from a JSON header - return NULL if not found
cJSON* json_hdr_try(cJSON *json, const char *field, int type, const char *path);
// Get values from a JSON header - die() if not found
cJSON* json_hdr_get(cJSON *json, const char *field, int type, const char *path);

long   json_hdr_demand_int( cJSON *root, const char *field, const char *path);
size_t json_hdr_demand_uint(cJSON *json, const char *field, const char *path);

#define json_hdr_get_graph(root,fpath) json_hdr_get(root,"graph",cJSON_Object,fpath)
#define json_hdr_get_paths(root,fpath) json_hdr_get(root,"paths",cJSON_Object,fpath)

size_t json_hdr_get_kmer_size(cJSON *root, const char *path);
size_t json_hdr_get_ncols(cJSON *json, const char *path);

// Get the number of non-ref samples in the graph
size_t json_hdr_get_nonref_ncols(cJSON *json, const char *path);
bool json_hdr_colour_is_ref(cJSON *json);

cJSON* json_hdr_get_curr_cmd(cJSON *json, const char *path);

#endif /* JSON_HDR_H_ */
