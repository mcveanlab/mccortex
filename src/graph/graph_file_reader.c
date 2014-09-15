#include "global.h"
#include "graph_file_reader.h"
#include "graph_format.h"
#include "db_node.h"
#include "cmd.h"
#include "file_util.h"

int graph_file_open(GraphFileReader *file, const char *path)
{
  return graph_file_open2(file, path, "r");
}

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open2(GraphFileReader *file, const char *input, const char *mode)
{
  GraphFileHeader *hdr = &file->hdr;
  FileFilter *fltr = &file->fltr;
  file_filter_open(fltr, input); // calls die() on error
  const char *path = fltr->path.b;

  // Stat will fail on streams, so file_size and num_of_kmers with both be -1
  struct stat st;
  file->file_size = -1;
  file->num_of_kmers = -1;

  if(strcmp(input,"-") != 0) {
    if(stat(path, &st) == 0) file->file_size = st.st_size;
    else warn("Couldn't get file size: %s", futil_outpath_str(path));
  }

  file->fh = futil_fopen(path, mode);
  file->hdr_size = graph_file_read_header(file->fh, hdr, path);

  file_filter_set_cols(fltr, hdr->num_of_cols);

  // Check we can handle the kmer size
  db_graph_check_kmer_size(file->hdr.kmer_size, file->fltr.path.b);

  size_t bytes_per_kmer, bytes_remaining;

  // If reading from STDIN we don't know file size
  if(file->file_size != -1)
  {
    // File header checks
    // Get number of kmers
    bytes_per_kmer = sizeof(BinaryKmer) +
                     hdr->num_of_cols * (sizeof(Covg) + sizeof(Edges));
    bytes_remaining = (size_t)(file->file_size - file->hdr_size);
    file->num_of_kmers = (bytes_remaining / bytes_per_kmer);

    if(bytes_remaining % bytes_per_kmer != 0) {
      warn("Truncated graph file: %s [bytes per kmer: %zu "
           "remaining: %zu; fsize: %zu; header: %zu; nkmers: %zu]",
           path, bytes_per_kmer, bytes_remaining,
           (size_t)file->file_size, (size_t)file->hdr_size,
           (size_t)file->num_of_kmers);
    }
  }

  return 1;
}

// Close file
void graph_file_close(GraphFileReader *file)
{
  if(file->fh) fclose(file->fh);
  file_filter_close(&file->fltr);
  graph_header_dealloc(&file->hdr);
  memset(file, 0, sizeof(*file));
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// be sure to zero covgs, edges before reading in
bool graph_file_read(const GraphFileReader *file,
                     BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  // status("Header colours: %u", file->hdr.num_of_cols);
  Covg kmercovgs[file->hdr.num_of_cols];
  Edges kmeredges[file->hdr.num_of_cols];
  size_t i, from, into;
  const FileFilter *fltr = &file->fltr;

  if(!graph_file_read_kmer(file->fh, &file->hdr, fltr->path.b,
                           bkmer, kmercovgs, kmeredges)) return false;

  for(i = 0; i < fltr->ncols; i++) {
    from = file_filter_fromcol(fltr, i);
    into = file_filter_intocol(fltr, i);
    covgs[into] = SAFE_ADD_COVG(covgs[into], kmercovgs[from]);
    edges[into] |= kmeredges[from];
  }

  return true;
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
// @ncols is file_filter_into_ncols(&file->fltr)
bool graph_file_read_reset(const GraphFileReader *file, size_t ncols,
                           BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  memset(covgs, 0, ncols*sizeof(Covg));
  memset(edges, 0, ncols*sizeof(Edges));
  return graph_file_read(file, bkmer, covgs, edges);
}

// Returns true if one or more files passed loads data into colour
bool graph_file_is_colour_loaded(size_t colour, const GraphFileReader *files,
                                 size_t num_files)
{
  size_t i;
  for(i = 0; i < num_files; i++) {
    if(file_filter_iscolloaded(&files[i].fltr, colour))
      return true;
  }
  return false;
}

// if one of the files is reading from stdin, sum_kmers_ptr is set to 0
// `max_cols_ptr` is used to return the most colours being loaded from a single file
// returns the number of colours being loaded in total
size_t graph_files_open(const char **graph_paths,
                        GraphFileReader *gfiles, size_t num_gfiles,
                        size_t *max_kmers_ptr, size_t *sum_kmers_ptr)
{
  size_t i, ctx_max_kmers = 0, ctx_sum_kmers = 0;
  bool ctx_uses_stdin = false;
  size_t ncols = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    memset(&gfiles[i], 0, sizeof(GraphFileReader));
    graph_file_open(&gfiles[i], graph_paths[i]);

    if(gfiles[0].hdr.kmer_size != gfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                      gfiles[0].hdr.kmer_size, gfiles[i].hdr.kmer_size);
    }

    file_filter_shift_cols(&gfiles[i].fltr, ncols);
    ncols = MAX2(ncols, file_filter_into_ncols(&gfiles[i].fltr));

    ctx_max_kmers = MAX2(ctx_max_kmers, graph_file_nkmers(&gfiles[i]));
    ctx_sum_kmers += graph_file_nkmers(&gfiles[i]);
    ctx_uses_stdin |= file_filter_isstdin(&gfiles[i].fltr);
  }

  if(ctx_uses_stdin) ctx_sum_kmers = SIZE_MAX;

  *max_kmers_ptr = ctx_max_kmers;
  *sum_kmers_ptr = ctx_sum_kmers;

  return ncols;
}
