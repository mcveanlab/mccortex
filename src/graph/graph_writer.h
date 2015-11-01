#ifndef GRAPH_WRITER_H_
#define GRAPH_WRITER_H_

#include <inttypes.h>
#include "db_graph.h"
#include "graph_format.h"
#include "graph_file_reader.h"

//
// Write graphs files to disk and merge graphs files on disk
//

/*!
  Write kmers from the graph to a file. The file header should already have been
  written.
  @return Number of bytes written
 */
size_t graph_write_empty(const dBGraph *db_graph, FILE *fh, size_t num_of_cols);

/*!
  Overwrite kmers in an existing file.
  @param first_graphcol first colour in the dBGraph to read from
  @param first_filecol first colour in the file to write into
  @param ngraphcols Number of colours to write to file
  @param nfilecols Total number of colours in file
  @param mmap_ptr Memory mapped file pointer
  @param hdrsize Size of file header i.e. byte pos of first kmer in file
 */
void graph_writer_update_mmap_kmers(const dBGraph *db_graph,
                                    size_t first_graphcol, size_t ngraphcols,
                                    size_t first_filecol, size_t nfilecols,
                                    char *mmap_ptr, size_t hdrsize);

/*!
  Overwrite kmers in an existing file.
  @param first_graphcol  first colour in the dBGraph to read from
  @param first_filecol   first colour in the file to write into
  @param ngraphcols      number of colours to write to file
  @param nfilecols       total number of colours in file
  @param hdrsize         file header size i.e. byte pos of first kmer in file
  @param fh              file handle opened that we can fwrite/fseek
  @param path            path to the file (for error messages)
 */
void graph_writer_update_file_kmers(const dBGraph *db_graph,
                                    size_t first_graphcol, size_t ngraphcols,
                                    size_t first_filecol, size_t nfilecols,
                                    size_t hdrsize, FILE *fh, const char *path);

// Returns number of bytes written
size_t graph_write_header(FILE *fh, const GraphFileHeader *header);

size_t graph_write_kmer(FILE *fh, size_t num_bkmer_words, size_t num_cols,
                        const BinaryKmer bkmer, const Covg *covgs,
                        const Edges *edges);

// Dump all kmers with all colours to given file. Returns num of kmers written
size_t graph_write_all_kmers(FILE *fh, const dBGraph *db_graph);

// If you don't want to/care about graph_info, pass in NULL
// If you want to print all nodes pass condition as NULL
// start_col is ignored unless colours is NULL
// returns number of nodes dumped
uint64_t graph_writer_save_mkhdr(const char *path, const dBGraph *graph,
                                 uint32_t version,
                                 const Colour *colours, Colour start_col,
                                 size_t num_of_cols);

// Pass your own header
uint64_t graph_writer_save(const char *path, const dBGraph *db_graph,
                           const GraphFileHeader *header, size_t intocol,
                           const Colour *colours, Colour start_col,
                           size_t num_of_cols);

void graph_writer_print_status(uint64_t nkmers, size_t ncols,
                               const char *path, uint32_t version);

//
// Merging, filtering, combining graph files
//

// Load a kmer and write to a file one kmer at a time
// Optionally filter a against the graph currently loaded
//   (i.e. only keep nodes and edges that are in the graph)
// Same functionality as graph_writer_merge, but faster if dealing with only one
// input file. Reads in and dumps one kmer at a time
size_t graph_writer_stream(const char *out_ctx_path, GraphFileReader *file,
                           const dBGraph *db_graph, const GraphFileHeader *hdr,
                           const Edges *only_load_if_in_edges);

size_t graph_writer_stream_mkhdr(const char *out_ctx_path, GraphFileReader *file,
                                 const dBGraph *db_graph,
                                 const Edges *only_load_if_in_edges,
                                 const char *intersect_gname);

size_t graph_writer_merge(const char *out_ctx_path,
                          GraphFileReader *files, size_t num_files,
                          bool kmers_loaded, bool colours_loaded,
                          const Edges *only_load_if_in_edges,
                          GraphFileHeader *hdr, dBGraph *db_graph);

// if intersect only load kmers that are already in the hash table
// returns number of kmers written
size_t graph_writer_merge_mkhdr(const char *out_ctx_path,
                                GraphFileReader *files, size_t num_files,
                                bool kmers_loaded, bool colours_loaded,
                                const Edges *only_load_if_in_edges,
                                const char *intersect_gname, dBGraph *db_graph);

#endif /* GRAPH_WRITER_H_ */
