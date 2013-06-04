#ifndef ADD_READ_PATH_H_
#define ADD_READ_PATH_H_

#include "file_reader.h"

void add_read_paths_to_graph(const char *se_list,
                             const char *pe_list1, const char *pe_list2,
                             Colour seq_colour,
                             const char *colour_list,
                             Colour col_list_first_colour,
                             SeqLoadingPrefs prefs);

#endif /* ADD_READ_PATH_H_ */
