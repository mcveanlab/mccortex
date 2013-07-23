#ifndef BUBBLE_CALLER_H_
#define BUBBLE_CALLER_H_

#include "db_graph.h"

void invoke_bubble_caller(const dBGraph *db_graph, const char* out_file,
                          int num_threads, char **tmp_paths);

#endif /* BUBBLE_CALLER_H_ */
