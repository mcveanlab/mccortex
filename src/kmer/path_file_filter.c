#include "global.h"
#include "path_file_filter.h"
#include "path_format.h"
#include "util.h"

const PathFileHeader INIT_PATH_FILE_HDR = INIT_PATH_FILE_HDR_MACRO;
const PathFileReader INIT_PATH_READER = INIT_PATH_READER_MACRO;

int path_file_open(PathFileReader *file, char *path, bool fatal)
{
  return path_file_open2(file, path, fatal, "r");
}

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new PathFileReader and returns 1
int path_file_open2(PathFileReader *file, char *path, bool fatal,
                    const char *mode)
{
  PathFileHeader *hdr = &file->hdr;
  FileFilter *fltr = &file->fltr;

  if(!file_filter_alloc(fltr, path, mode, fatal)) return 0;
  setvbuf(fltr->fh, NULL, _IOFBF, CTP_BUF_SIZE);

  file->hdr_size = paths_file_read_header(fltr->fh, hdr, fatal, fltr->file_path.buff);
  if(file->hdr_size == -1) return -1;

  file_filter_set_cols(fltr, hdr->num_of_cols);

  // Check we can handle the kmer size
  db_graph_check_kmer_size(file->hdr.kmer_size, file->fltr.file_path.buff);

  // Check file length
  size_t file_len = (size_t)file->hdr_size + (size_t)hdr->num_path_bytes +
                    hdr->num_kmers_with_paths *
                      (NUM_BKMER_WORDS*sizeof(uint64_t) + sizeof(uint64_t));

  if(file_len != (size_t)file->fltr.file_size) {
    warn("Corrupted file? Sizes don't match up "
         "[hdr:%zu exp:%zu actual:%zu path: %s]",
         (size_t)file->hdr_size, file_len,
         (size_t)file->fltr.file_size, file->fltr.file_path.buff);
  }

  return 1;
}

// File header checks
void path_file_load_check(const PathFileReader *file, const dBGraph *db_graph)
{
  const FileFilter *fltr = &file->fltr;
  const PathFileHeader *hdr = &file->hdr;

  // Conservative set of tests to see if we can hold the data from a path file
  if(file->hdr.kmer_size != db_graph->kmer_size) {
    die("Kmer sizes do not match between graph and path file [%s]",
        fltr->orig_path.buff);
  }

  if(hdr->num_path_bytes > db_graph->pdata.size) {
    char mem_str[100]; bytes_to_str(hdr->num_path_bytes, 1, mem_str);
    die("Not enough memory allocated to store paths [mem: %s path: %s]",
        mem_str, fltr->orig_path.buff);
  }

  if(db_graph->ht.num_kmers > 0 &&
     db_graph->ht.num_kmers < hdr->num_kmers_with_paths)
  {
    warn("Graph has fewer kmers than paths file");
  }

  // More checks to ensure we can load and use a path file along with a graph file
  if(path_file_usedcols(file) > db_graph->num_of_cols) {
    die("Number of colours in path file is greater than in the graph %zu > %zu [%s]",
        path_file_usedcols(file), db_graph->num_of_cols, fltr->orig_path.buff);
  }

  // Check sample names match
  size_t i, intocol, fromcol;
  for(i = 0; i < hdr->num_of_cols; i++)
  {
    intocol = file_filter_intocol(fltr, i);
    fromcol = file_filter_fromcol(fltr, i);
    char *gname = db_graph->ginfo[intocol].sample_name.buff;
    char *pname = hdr->sample_names[fromcol].buff;

    if(strcmp(pname, "noname") != 0 && strcmp(gname, pname) != 0) {
      die("Graph/path sample names do not match [%zu->%zu] '%s' vs '%s'",
          i, intocol, gname, pname);
    }
  }
}

void path_file_set_graph_sample_names(const PathFileReader *file,
                                      dBGraph *db_graph)
{
  const FileFilter *fltr = &file->fltr;
  const PathFileHeader *hdr = &file->hdr;

  size_t i, intocol, fromcol;
  StrBuf *gname, *pname;
  for(i = 0; i < hdr->num_of_cols; i++)
  {
    intocol = file_filter_intocol(fltr, i);
    fromcol = file_filter_fromcol(fltr, i);
    gname = &db_graph->ginfo[intocol].sample_name;
    pname = &hdr->sample_names[fromcol];

    if(strcmp(gname->buff, "undefined") == 0)
      strbuf_set(gname, pname->buff);
    else if(strcmp(gname->buff, pname->buff) != 0) {
      die("Graph/path sample names do not match [%zu->%zu] '%s' vs '%s'",
          fromcol, intocol, gname->buff, pname->buff);
    }
  }
}

void path_file_set_header_sample_names(const PathFileReader *file,
                                       PathFileHeader *hdr1)
{
  const FileFilter *fltr = &file->fltr;
  const PathFileHeader *hdr0 = &file->hdr;

  size_t i, intocol, fromcol;
  const StrBuf *name0; StrBuf *name1;

  for(i = 0; i < fltr->ncols; i++)
  {
    intocol = file_filter_intocol(fltr, i);
    fromcol = file_filter_fromcol(fltr, i);
    name1 = &hdr1->sample_names[intocol];
    name0 = &hdr0->sample_names[fromcol];
    strbuf_set(name1, name0->buff);
  }
}

// Close file
void path_file_close(PathFileReader *file)
{
  file_filter_close(&file->fltr);
}

// calls file_filter_dealloc which will close file if needed
void path_file_dealloc(PathFileReader *file)
{
  file_filter_dealloc(&file->fltr);
  paths_header_dealloc(&file->hdr);
}
