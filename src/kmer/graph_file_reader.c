#include "global.h"
#include "graph_file_reader.h"
#include "graph_format.h"
#include "cmd.h"

const GraphFileHeader INIT_GRAPH_FILE_HDR = INIT_GRAPH_FILE_HDR_MACRO;
const GraphFileReader INIT_GRAPH_READER = INIT_GRAPH_READER_MACRO;

int graph_file_open(GraphFileReader *file, char *path, bool fatal)
{
  return graph_file_open2(file, path, fatal, "r");
}

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open2(GraphFileReader *file, char *path, bool fatal,
                     const char *mode)
{
  GraphFileHeader *hdr = &file->hdr;
  FileFilter *fltr = &file->fltr;

  if(!file_filter_open(fltr, path, mode, fatal)) return 0;
  setvbuf(fltr->fh, NULL, _IOFBF, CTX_BUF_SIZE);

  file->hdr_size = graph_file_read_header(fltr->fh, hdr, fatal, fltr->file_path.buff);
  if(file->hdr_size == -1) return -1;

  file_filter_set_cols(fltr, hdr->num_of_cols);

  // Check we can handle the kmer size
  db_graph_check_kmer_size(file->hdr.kmer_size, file->fltr.file_path.buff);

  size_t bytes_per_kmer, bytes_remaining, nkmers = 0;

  // If reading from STDIN we don't know file size
  if(fltr->file_size != -1)
  {
    // File header checks
    // Get number of kmers
    bytes_per_kmer = sizeof(BinaryKmer) +
                     hdr->num_of_cols * (sizeof(Covg) + sizeof(Edges));
    bytes_remaining = (size_t)(fltr->file_size - file->hdr_size);
    nkmers = (bytes_remaining / bytes_per_kmer);

    if(bytes_remaining % bytes_per_kmer != 0) {
      warn("Truncated graph file: %s [bytes per kmer: %zu "
           "remaining: %zu; fsize: %zu; header: %zu; nkmers: %zu]",
           fltr->file_path.buff, bytes_per_kmer, bytes_remaining,
           (size_t)fltr->file_size, (size_t)file->hdr_size, nkmers);
    }
  }

  file->num_of_kmers = nkmers;

  return 1;
}

// Close file
void graph_file_close(GraphFileReader *file)
{
  file_filter_close(&file->fltr);
  graph_header_dealloc(&file->hdr);
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
bool graph_file_read(const GraphFileReader *file,
                     BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  // status("Header colours: %u", file->hdr.num_of_cols);
  Covg kmercovgs[file->hdr.num_of_cols];
  Edges kmeredges[file->hdr.num_of_cols];
  size_t i;
  const FileFilter *fltr = &file->fltr;

  if(!graph_file_read_kmer(fltr->fh, &file->hdr, fltr->file_path.buff,
                           bkmer, kmercovgs, kmeredges)) return false;

  // covgs += file->intocol;
  // edges += file->intocol;

  if(fltr->flatten) {
    covgs[0] = 0;
    edges[0] = 0;
    for(i = 0; i < fltr->ncols; i++) {
      covgs[0] += kmercovgs[fltr->cols[i]];
      edges[0] |= kmeredges[fltr->cols[i]];
    }
  }
  else {
    for(i = 0; i < fltr->ncols; i++) {
      covgs[i] = kmercovgs[fltr->cols[i]];
      edges[i] = kmeredges[fltr->cols[i]];
    }
  }

  return true;
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
size_t graph_files_open(char **graph_paths,
                        GraphFileReader *gfiles, size_t num_gfiles,
                        size_t *max_kmers_ptr, size_t *sum_kmers_ptr)
{
  size_t i, ctx_max_kmers = 0, ctx_sum_kmers = 0;
  bool ctx_uses_stdin = false;
  size_t ncols = 0;

  for(i = 0; i < num_gfiles; i++)
  {
    gfiles[i] = INIT_GRAPH_READER;
    graph_file_open(&gfiles[i], graph_paths[i], true);

    if(gfiles[0].hdr.kmer_size != gfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                      gfiles[0].hdr.kmer_size, gfiles[i].hdr.kmer_size);
    }

    file_filter_update_intocol(&gfiles[i].fltr, ncols);
    ncols = MAX2(ncols, graph_file_usedcols(&gfiles[i]));

    ctx_max_kmers = MAX2(ctx_max_kmers, gfiles[i].num_of_kmers);
    ctx_sum_kmers += gfiles[i].num_of_kmers;
    ctx_uses_stdin |= file_filter_isstdin(&gfiles[i].fltr);
  }

  if(ctx_uses_stdin) ctx_sum_kmers = SIZE_MAX;

  *max_kmers_ptr = ctx_max_kmers;
  *sum_kmers_ptr = ctx_sum_kmers;

  return ncols;
}
