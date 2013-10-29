#ifndef GRAPH_FILE_FILTER_H_
#define GRAPH_FILE_FILTER_H_

#include "graph_typedef.h"

typedef struct
{
  // uint32_t file_ncols;
  uint32_t version, kmer_size, num_of_bitfields, num_of_cols, max_col;
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

#define INIT_GRAPH_FILE_HDR {.capacity = 0}

const GraphFileReader INIT_GRAPH_READER;

// 4MB buffer
#define CTX_BUF_SIZE (4UL<<20)

#define graph_file_outncols(gfr) ((gfr)->flatten ? 1 : (gfr)->ncols)

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open(GraphFileReader *file, char *path, boolean fatal);

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

// Print file filter description
void graph_file_status(const GraphFileReader *file);

#endif /* GRAPH_FILE_FILTER_H_ */
