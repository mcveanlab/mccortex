#ifndef GRAPH_FILE_FILTER_H_
#define GRAPH_FILE_FILTER_H_

#include "graph_typedef.h"

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
  off_t file_size, hdr_size;
  StrBuf path; // file path without modifiers e.g. "2:in.ctx:0,3" path == "in.ctx"
  FILE *fh;
  uint32_t intocol, ncols, *cols, ncolscap;
  boolean flatten; // Merge all colours into intocol
} GraphFileReader;

#define INIT_GRAPH_FILE_HDR_MACRO {                    \
  .version = CTX_GRAPH_FILEFORMAT,                     \
  .kmer_size = 0, .num_of_bitfields = NUM_BKMER_WORDS, \
  .num_of_cols = 0, .num_of_kmers = 0, .capacity = 0}

#define INIT_GRAPH_READER_MACRO {                            \
  .hdr = INIT_GRAPH_FILE_HDR_MACRO, .file_size = 0, .hdr_size = 0, \
  .path = {.buff = NULL}, .fh = NULL,                        \
  .intocol = 0, .ncols = 0, .cols = NULL, .ncolscap = 0, .flatten = false}

const GraphFileHeader INIT_GRAPH_FILE_HDR;
const GraphFileReader INIT_GRAPH_READER;

// 4MB buffer
#define CTX_BUF_SIZE (4UL<<20)

#define graph_file_outncols(reader) ((reader)->flatten ? 1 : (reader)->ncols)
#define graph_file_intocol(reader,col) ((reader)->intocol + (!(reader)->flatten)*(col))

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

// Return true if all colours are being loaded once in their original order
boolean graph_file_no_filter(const GraphFileReader *file);

// Print file filter description
void graph_file_status(const GraphFileReader *file);

#endif /* GRAPH_FILE_FILTER_H_ */
