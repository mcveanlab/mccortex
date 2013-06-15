#ifndef SHADED_CALLER_H_
#define SHADED_CALLER_H_

#include "db_graph.h"

void invoke_shaded_bubble_caller(const dBGraph *db_graph, const char* out_file,
                                 int num_threads, char **tmp_paths);

#endif /* SHADED_CALLER_H_ */
