#ifndef GRAPH_WRITER_H_
#define GRAPH_WRITER_H_

#include <inttypes.h>
#include "db_graph.h"
#include "graph_format.h"
#include "graph_file_reader.h"

//
// Write graphs files to disk and merge graphs files on disk
//

// Construct graph header
// Free with graph_header_free(hdr)
GraphFileHeader* graph_writer_mkhdr(const dBGraph *db_graph,
                                    const FileFilter *fltr,
                                    size_t filencols);

// Returns number of bytes written
size_t graph_write_header(FILE *fh, const GraphFileHeader *hdr);

// Write a single kmer with its edges and coverages for colours 0..ncols-1 to
// 0..ncols-1 in the file. `covgs` and `edges` should be of length `filencols`.
// Returns number of bytes written
size_t graph_write_kmer(FILE *fh, size_t filencols,
                        const BinaryKmer bkmer,
                        const Covg *covgs,
                        const Edges *edges);

// Dump all kmers with all colours to given file.
// Do not use filters to re-arrange colours
// `sort_kmer` if true, sort kmers before writing. Uses extra memory:
//   requires sizeof(hkey_t) * num_kmers memory which it allocates and frees
// Returns num of kmers written
size_t graph_write_all_kmers_direct(FILE *fh, const dBGraph *db_graph,
                                    bool sort_kmers, const GraphFileHeader *hdr);

// Dump all kmers with all colours to given file.
// Filter kmers and re-arrange colours
// `sort_kmer` if true, sort kmers before writing. Uses extra memory:
//   requires sizeof(hkey_t) * num_kmers memory which it allocates and frees
// Returns num of kmers written
size_t graph_write_all_kmers_filtered(FILE *fh, const dBGraph *db_graph,
                                      bool sort_kmers, const GraphFileHeader *hdr,
                                      const FileFilter *fltr);

// Pass your own header, ncols in file taken from the header
// If sort_kmers is true, save kmers in lexigraphical order
// returns number of nodes written out
uint64_t graph_writer_save(const char *path, const dBGraph *db_graph,
                           const GraphFileHeader *hdr, bool sort_kmers,
                           const FileFilter *fltr);

// filencols must be <= db_graph->num_of_cols
// returns number of nodes dumped
uint64_t graph_writer_save_mkhdr(const char *path, const dBGraph *db_graph,
                                 bool sort_kmers, size_t filencols);

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
                          GraphFileHeader *hdr, bool sort_kmers,
                          dBGraph *db_graph);

// if intersect only load kmers that are already in the hash table
// returns number of kmers written
size_t graph_writer_merge_mkhdr(const char *out_ctx_path,
                                GraphFileReader *files, size_t num_files,
                                bool kmers_loaded, bool colours_loaded,
                                const Edges *only_load_if_in_edges,
                                const char *intersect_gname,
                                bool sort_kmers, dBGraph *db_graph);

#endif /* GRAPH_WRITER_H_ */
