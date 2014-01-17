#ifndef GRAPH_FILE_FILTER_H_
#define GRAPH_FILE_FILTER_H_

#include "graph_typedef.h"
#include "file_filter.h"
#include "binary_kmer.h"

// Graph (.ctx)
typedef struct
{
  uint32_t version, kmer_size, num_of_bitfields, num_of_cols;
  uint64_t num_of_kmers;
  GraphInfo *ginfo; // Cleaning info etc for each colour
  size_t capacity; // number of ginfo objects malloc'd
} GraphFileHeader;

typedef struct
{
  GraphFileHeader hdr;
  off_t hdr_size;
  FileFilter fltr; // colour filter
} GraphFileReader;

#define INIT_GRAPH_FILE_HDR_MACRO {                    \
  .version = CTX_GRAPH_FILEFORMAT,                     \
  .kmer_size = 0, .num_of_bitfields = NUM_BKMER_WORDS, \
  .num_of_cols = 0, .num_of_kmers = 0, .capacity = 0}

#define INIT_GRAPH_READER_MACRO {                  \
  .hdr = INIT_GRAPH_FILE_HDR_MACRO, .hdr_size = 0, \
  .fltr = INIT_FILE_FILTER_MACRO}

const GraphFileHeader INIT_GRAPH_FILE_HDR;
const GraphFileReader INIT_GRAPH_READER;

// 4MB buffer
#define CTX_BUF_SIZE (4UL<<20)

// #define graph_file_outncols(reader) file_filter_outncols(reader)
// #define graph_file_intocol(reader,col) file_filter_intocol(reader,col)

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
int graph_file_open(GraphFileReader *file, char *path, boolean fatal);

// mode is "r", "r+" etc.
int graph_file_open2(GraphFileReader *file, char *path, boolean fatal,
                     const char *mode);

// Close file
void graph_file_close(GraphFileReader *file);

// Release all memory (also calls close)
void graph_file_dealloc(GraphFileReader *file);

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// Beware: this function does not use file.intocol so you may wish to pass:
//    graph_file_read(file, &bkmer, covgs+file.intocol, edges+file.intocol);
boolean graph_file_read(const GraphFileReader *file,
                        BinaryKmer *bkmer, Covg *covgs, Edges *edges);

// Returns true if one or more files passed loads data into colour
boolean graph_file_is_colour_loaded(size_t colour, const GraphFileReader *files,
                                    size_t num_files);

#endif /* GRAPH_FILE_FILTER_H_ */
