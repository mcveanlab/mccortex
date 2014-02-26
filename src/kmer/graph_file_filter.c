#include "global.h"
#include "graph_file_filter.h"
#include "graph_format.h"

const GraphFileHeader INIT_GRAPH_FILE_HDR = INIT_GRAPH_FILE_HDR_MACRO;
const GraphFileReader INIT_GRAPH_READER = INIT_GRAPH_READER_MACRO;

int graph_file_open(GraphFileReader *file, char *path, boolean fatal)
{
  return graph_file_open2(file, path, fatal, "r");
}

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GraphFileReader and returns 1
int graph_file_open2(GraphFileReader *file, char *path, boolean fatal,
                     const char *mode)
{
  GraphFileHeader *hdr = &file->hdr;
  FileFilter *fltr = &file->fltr;

  if(!file_filter_alloc(fltr, path, mode, fatal)) return 0;
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

    if(hdr->version > 6 && hdr->num_of_kmers != nkmers) {
      warn("File size and number of kmers do not match: %s [bytes per kmer: %zu "
           "remaining: %zu; fsize: %zu; header: %zu; expect: %zu; got: %zu]",
           fltr->file_path.buff, bytes_per_kmer, bytes_remaining,
           (size_t)fltr->file_size, (size_t)file->hdr_size,
           (size_t)hdr->num_of_kmers, nkmers);
    }

    if(bytes_remaining % bytes_per_kmer != 0) {
      warn("Truncated graph file: %s [bytes per kmer: %zu "
           "remaining: %zu; fsize: %zu; header: %zu; nkmers: %zu]",
           fltr->file_path.buff, bytes_per_kmer, bytes_remaining,
           (size_t)fltr->file_size, (size_t)file->hdr_size, nkmers);
    }
  }

  hdr->num_of_kmers = nkmers;

  return 1;
}

// Close file
void graph_file_close(GraphFileReader *file)
{
  file_filter_close(&file->fltr);
}

// calls file_filter_dealloc which will close file if needed
void graph_file_dealloc(GraphFileReader *file)
{
  file_filter_dealloc(&file->fltr);
  graph_header_dealloc(&file->hdr);
}

// Read a kmer from the file
// returns true on success, false otherwise
// prints warnings if dirty kmers in file
boolean graph_file_read(const GraphFileReader *file,
                        BinaryKmer *bkmer, Covg *covgs, Edges *edges)
{
  // status("Header colours: %u", file->hdr.num_of_cols);
  Covg kmercovgs[file->hdr.num_of_cols];
  Edges kmeredges[file->hdr.num_of_cols];
  size_t i;
  const FileFilter *fltr = &file->fltr;

  if(!graph_file_read_kmer(fltr->fh, &file->hdr, fltr->file_path.buff,
                           bkmer, kmercovgs, kmeredges)) return 0;

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

  return 1;
}

// Returns true if one or more files passed loads data into colour
boolean graph_file_is_colour_loaded(size_t colour, const GraphFileReader *files,
                                    size_t num_files)
{
  size_t i;
  for(i = 0; i < num_files; i++) {
    if(file_filter_iscolloaded(&files[i].fltr, colour))
      return true;
  }
  return false;
}
