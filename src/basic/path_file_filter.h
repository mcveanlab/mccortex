#ifndef PATH_FILE_FILTER_H_
#define PATH_FILE_FILTER_H_

#include "graph_typedef.h"
#include "file_filter.h"
#include "graph_format.h"

typedef struct
{
  PathFileHeader hdr;
  off_t hdr_size;
  FileFilter fltr; // colour filter
} PathFileReader;

#define INIT_PATH_FILE_HDR_MACRO {                \
  .version = CTX_PATH_FILEFORMAT,                 \
  .kmer_size = 0, .num_of_paths = 0,              \
  .num_path_bytes = 0, .num_kmers_with_paths = 0, \
  .sample_names = NULL, .capacity = 0}

#define INIT_PATH_READER_MACRO {                  \
  .hdr = INIT_PATH_FILE_HDR_MACRO, .hdr_size = 0, \
  .fltr = INIT_FILE_FILTER_MACRO}

const PathFileHeader INIT_PATH_FILE_HDR;
const PathFileReader INIT_PATH_READER;

// 4MB buffer
#define CTP_BUF_SIZE (4UL<<20)

#define path_file_outncols(rdr) file_filter_outncols(&(rdr)->fltr)
#define path_file_intocol(rdr,col) file_filter_intocol(&(rdr)->fltr,col)
#define path_file_usedcols(rdr) file_filter_usedcols(&(rdr)->fltr)

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new PathFileReader and returns 1
int path_file_open(PathFileReader *file, char *path, boolean fatal);

// mode is "r", "r+" etc.
int path_file_open2(PathFileReader *file, char *path, boolean fatal,
                    const char *mode);

// File header checks
void path_file_load_check(const PathFileReader *file, const dBGraph *db_graph);

// Close file
void path_file_close(PathFileReader *file);

// Release all memory (also calls close)
void path_file_dealloc(PathFileReader *file);

#endif /* PATH_FILE_FILTER_H_ */
