#ifndef GRAPH_FILE_READER_H_
#define GRAPH_FILE_READER_H_

#include "cortex_types.h"
#include "file_filter.h"
#include "binary_kmer.h"
#include "graph_info.h"

// Graph (.ctx)
typedef struct
{
  uint32_t version, kmer_size, num_of_bitfields, num_of_cols;
  GraphInfo *ginfo; // Cleaning info etc for each colour
  size_t capacity; // number of ginfo objects malloc'd
} GraphFileHeader;

typedef struct
{
  FileFilter fltr;
  GraphFileHeader hdr;
  off_t hdr_size;
  uint64_t num_of_kmers; // only set if reading from file (i.e. not stream)
} GraphFileReader;

#define INIT_GRAPH_FILE_HDR_MACRO {                    \
  .version = CTX_GRAPH_FILEFORMAT,                     \
  .kmer_size = 0, .num_of_bitfields = NUM_BKMER_WORDS, \
  .num_of_cols = 0, .ginfo = NULL, .capacity = 0}

#define INIT_GRAPH_READER_MACRO {                   \
  .fltr = INIT_FILE_FILTER_MACRO, .num_of_kmers = 0,\
  .hdr = INIT_GRAPH_FILE_HDR_MACRO, .hdr_size = 0}

const GraphFileHeader INIT_GRAPH_FILE_HDR;
const GraphFileReader INIT_GRAPH_READER;

#include "objbuf_macro.h"
create_objbuf(gfile_buf, GraphFileBuffer, GraphFileReader);

#define graph_file_outncols(rdr) file_filter_outncols(&(rdr)->fltr)
#define graph_file_intocol(rdr,col) file_filter_intocol(&(rdr)->fltr,col)
#define graph_file_fromcol(rdr,col) file_filter_fromcol(&(rdr)->fltr,col)
#define graph_file_usedcols(rdr) file_filter_usedcols(&(rdr)->fltr)

// File path plus colour specification e.g. in.ctx:1,3
#define graph_file_orig_path(rdr) ((rdr)->fltr.orig_path.buff)
// Just the file path e.g. in.ctx
#define graph_file_file_path(rdr) ((rdr)->fltr.orig_path.buff)

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open(GraphFileReader *file, char *path, bool fatal);

// mode is "r", "r+" etc.
int graph_file_open2(GraphFileReader *file, char *path, bool fatal,
                     const char *mode);

// Close file, release all memory
void graph_file_close(GraphFileReader *file);

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// Beware: this function does not use file.intocol so you may wish to pass:
//    graph_file_read(file, &bkmer, covgs+file.intocol, edges+file.intocol);
bool graph_file_read(const GraphFileReader *file,
                        BinaryKmer *bkmer, Covg *covgs, Edges *edges);

// Returns true if one or more files passed loads data into colour
bool graph_file_is_colour_loaded(size_t colour, const GraphFileReader *files,
                                    size_t num_files);

// if one of the files is reading from stdin, sum_kmers_ptr is set to 0
// `max_cols_ptr` is used to return the most colours being loaded from a single file
// returns the number of colours being loaded in total
size_t graph_files_open(char **graph_paths,
                        GraphFileReader *gfiles, size_t num_gfiles,
                        size_t *max_kmers_ptr, size_t *sum_kmers_ptr);

#endif /* GRAPH_FILE_READER_H_ */
