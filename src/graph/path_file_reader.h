#ifndef PATH_FILE_READER_H_
#define PATH_FILE_READER_H_

#include "cortex_types.h"
#include "file_filter.h"
#include "graph_format.h"

// Path (.ctp)
typedef struct
{
  uint32_t version, kmer_size, num_of_cols;
  uint64_t num_of_paths, num_path_bytes, num_kmers_with_paths;
  StrBuf *sample_names;
  size_t capacity; // how many sample_names have been malloc'd
} PathFileHeader;

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

#include "objbuf_macro.h"
create_objbuf(pfile_buf, PathFileBuffer, PathFileReader);

// 4MB buffer
#define CTP_BUF_SIZE (4UL<<20)

#define path_file_outncols(rdr) file_filter_outncols(&(rdr)->fltr)
#define path_file_intocol(rdr,col) file_filter_intocol(&(rdr)->fltr,col)
#define path_file_fromcol(rdr,col) file_filter_fromcol(&(rdr)->fltr,col)
#define path_file_usedcols(rdr) file_filter_usedcols(&(rdr)->fltr)

#define path_file_colset_bytes(rdr) (((rdr)->hdr.num_of_paths+3)/4)

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new PathFileReader and returns 1
int path_file_open(PathFileReader *file, char *path, bool fatal);

// mode is "r", "r+" etc.
int path_file_open2(PathFileReader *file, char *path, bool fatal,
                    const char *mode);

// Close file, release memory
void path_file_close(PathFileReader *file);

// File header checks
void path_file_load_check(const PathFileReader *file, const dBGraph *db_graph);

// Copy sample names form PathFileReader to graph
void path_file_set_graph_sample_names(const PathFileReader *file,
                                      dBGraph *db_graph);

void path_file_set_header_sample_names(const PathFileReader *file,
                                       PathFileHeader *hdr1);

#endif /* PATH_FILE_READER_H_ */
