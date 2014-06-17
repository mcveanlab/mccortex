#ifndef CALLER_OUTPUT_H_
#define CALLER_OUTPUT_H_

#include "db_graph.h"

void caller_gzprint_header(gzFile gzout, const char* out_file,
                           const char *format_str, const dBGraph *db_graph);

#endif /* CALLER_OUTPUT_H_ */
