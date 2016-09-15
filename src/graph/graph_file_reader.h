#ifndef GRAPH_FILE_READER_H_
#define GRAPH_FILE_READER_H_

// #include "cortex_types.h"
#include "graph_format.h"
#include "file_filter.h"
#include "binary_kmer.h"

//
// Read graph files from disk
//

typedef struct
{
  FILE *fh;
  StreamBuffer strm; // buffer input
  FileFilter fltr;
  GraphFileHeader hdr;
  off_t hdr_size, file_size;
  int64_t num_of_kmers; // set if reading from file (i.e. not stream) else -1
  bool error_zero_covg, error_missing_covg; // Whether we saw loading errors
} GraphFileReader;

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(gfile_buf, GraphFileBuffer, GraphFileReader);


#define graph_file_reset(rdr) memset(rdr, 0, sizeof(GraphFileReader))

// Returns 0 if not set instead of -1
#define graph_file_nkmers(rdr) ((uint64_t)MAX2((rdr)->num_of_kmers, 0))

// Get file offset of a given kmer
static inline off_t graph_file_offset(const GraphFileReader *gfr, size_t i)
{
  size_t s = sizeof(BinaryKmer)+gfr->fltr.srcncols*(sizeof(Covg)+sizeof(Edges));
  return gfr->hdr_size + s*i;
}

#define graph_file_is_buffered(file) ((file)->strm.b != NULL)
// Buffer size `bufsize` is in bytes
void graph_file_set_buffered(GraphFileReader *file, size_t bufsize);

int graph_file_fseek(GraphFileReader *file, off_t offset, int whence);
off_t graph_file_ftell(GraphFileReader *file);

// read `n` bytes from `file` into `ptr`
size_t graph_file_fread(GraphFileReader *file, void *ptr, size_t n);

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open(GraphFileReader *file, const char *path);

// mode is "r", "r+" etc.
// If there is no 'into' filter (e.g. 0,1:in.ctx 0,1 is 'into filter'),
// load into `into_offset..into_offset+N-1`
int graph_file_open2(GraphFileReader *file, const char *path, const char *mode,
                     bool usebuf, size_t into_offset);

// Close file, release all memory
void graph_file_close(GraphFileReader *file);

// Merge the header from this file into a new file
void graph_file_merge_header(GraphFileHeader *hdr, const GraphFileReader *file);

// Returns number of bytes read
size_t graph_file_read_raw(GraphFileReader *rdr,
                           BinaryKmer *bkmer, Covg *covgs, Edges *edges);

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// be sure to zero covgs, edges before reading in
bool graph_file_read(GraphFileReader *file,
                     BinaryKmer *bkmer, Covg *covgs, Edges *edges);

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
bool graph_file_read_reset(GraphFileReader *file,
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
