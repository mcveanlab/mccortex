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
  size_t capacity;
} GraphFileHeader;

typedef struct
{
  FILE *fh;
  FileFilter fltr;
  GraphFileHeader hdr;
  off_t hdr_size, file_size;
  int64_t num_of_kmers; // set if reading from file (i.e. not stream) else -1
} GraphFileReader;

#define graph_file_reset(rdr) memset(rdr, 0, sizeof(GraphFileReader))

// Returns 0 if not set instead of -1
#define graph_file_nkmers(rdr) ((uint64_t)MAX2((rdr)->num_of_kmers, 0))

#include "objbuf_macro.h"
create_objbuf(gfile_buf, GraphFileBuffer, GraphFileReader);

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open(GraphFileReader *file, const char *path);

// mode is "r", "r+" etc.
int graph_file_open2(GraphFileReader *file, const char *path, const char *mode);

// Close file, release all memory
void graph_file_close(GraphFileReader *file);

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// Beware: this function does not use file.intocol so you may wish to pass:
//    graph_file_read(file, &bkmer, covgs+file.intocol, edges+file.intocol);
bool graph_file_read(const GraphFileReader *file,
                     BinaryKmer *bkmer, Covg *covgs, Edges *edges);

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// @ncols is file_filter_into_ncols(&file->fltr)
bool graph_file_read_reset(const GraphFileReader *file, size_t ncols,
                           BinaryKmer *bkmer, Covg *covgs, Edges *edges);

// Returns true if one or more files passed loads data into colour
bool graph_file_is_colour_loaded(size_t colour, const GraphFileReader *files,
                                 size_t num_files);

// if one of the files is reading from stdin, sum_kmers_ptr is set to 0
// `max_cols_ptr` is used to return the most colours being loaded from a single file
// returns the number of colours being loaded in total
size_t graph_files_open(const char **graph_paths,
                        GraphFileReader *gfiles, size_t num_gfiles,
                        size_t *max_kmers_ptr, size_t *sum_kmers_ptr);

#endif /* GRAPH_FILE_READER_H_ */
