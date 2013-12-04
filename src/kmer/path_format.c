#include "global.h"
#include "path_format.h"
#include "db_graph.h"
#include "db_node.h"
#include "hash_table.h"
#include "path_store.h"
#include "util.h"
#include "file_util.h"

// Format:
// -- Header --
// "PATHS"<uint32_t:version><uint32_t:kmersize><uint32_t:num_of_cols>
// <uint64_t:num_of_paths><uint64_t:num_path_bytes><uint64_t:num_kmers_with_paths>
// -- Colours --
// <uint32_t:sname_len><uint8_t x sname_len:sample_name> x num_of_cols
// -- Data --
// <uint8_t:path_data>
// <binarykmer><uint64_t:path_index_fw><uint64_t:path_index_rv>

// out must be at least strlen(path)+4 bytes long
void paths_format_filename(const char *path, char *out)
{
  size_t len = strlen(path);
  memcpy(out, path, len);
  if(strcasecmp(path+len-4,".ctp") != 0)
  {
    if(strcasecmp(path+len-4,".ctx") != 0) len += 4;
    memcpy(out+len-4, ".ctp", 4);
  }
  out[len] = '\0';
}

void paths_header_alloc(PathFileHeader *h, size_t num_of_cols)
{
  size_t i, old_cap = h->capacity;

  if(h->capacity == 0) {
    h->sample_names = malloc2(num_of_cols * sizeof(StrBuf));
  }
  else if(num_of_cols > h->capacity) {
    h->sample_names = realloc2(h->sample_names, num_of_cols * sizeof(StrBuf));
  }

  for(i = old_cap; i < num_of_cols; i++) {
    strbuf_alloc(&h->sample_names[i], 256);
    strbuf_set(&h->sample_names[i], "noname");
  }

  h->capacity = MAX2(old_cap, num_of_cols);
}

void paths_header_dealloc(PathFileHeader *h)
{
  size_t i;
  if(h->capacity > 0) {
    for(i = 0; i < h->capacity; i++) strbuf_dealloc(&h->sample_names[i]);
    free(h->sample_names);
    h->capacity = 0;
  }
}

void paths_header_update(PathFileHeader *header, const PathStore *paths)
{
  header->num_of_paths = paths->num_of_paths;
  header->num_path_bytes = paths->next - paths->store;
  header->num_kmers_with_paths = paths->num_kmers_with_paths;
}

// Returns number of bytes read or -1 on error (if fatal is false)
int paths_file_read_header(FILE *fh, PathFileHeader *h,
                           boolean fatal, const char *path)
{
  int bytes_read = 0;
  char sig[6] = {0};

  SAFE_READ(fh, sig, 5, "PATHS", path, fatal);
  SAFE_READ(fh, &h->version, sizeof(uint32_t), "version", path, fatal);
  SAFE_READ(fh, &h->kmer_size, sizeof(uint32_t), "kmer_size", path, fatal);
  SAFE_READ(fh, &h->num_of_cols, sizeof(uint32_t), "num_of_cols", path, fatal);
  SAFE_READ(fh, &h->num_of_paths, sizeof(uint64_t), "num_of_paths", path, fatal);
  SAFE_READ(fh, &h->num_path_bytes, sizeof(uint64_t),
            "num_path_bytes", path, fatal);
  SAFE_READ(fh, &h->num_kmers_with_paths, sizeof(uint64_t),
            "num_kmers_with_paths", path, fatal);

  bytes_read += 5 + sizeof(uint32_t)*3 + sizeof(uint64_t)*3;

  // paths_header_alloc will only alloc or realloc only if it needs to
  paths_header_alloc(h, h->num_of_cols);

  // Read sample names
  size_t i;
  uint32_t len;
  StrBuf *sbuf;
  for(i = 0; i < h->num_of_cols; i++)
  {
    sbuf = h->sample_names + i;
    SAFE_READ(fh, &len, sizeof(uint32_t), "sample name length", path, fatal);
    strbuf_ensure_capacity(sbuf, len);
    SAFE_READ(fh, sbuf->buff, len, "sample name", path, fatal);
    sbuf->buff[sbuf->len = len] = '\0';
    bytes_read += sizeof(uint32_t) + len;
  }

  // Checks
  if(h->version < 1 || h->version > 1) {
    if(!fatal) return -1;
    die("file version not supported [version: %u; path: %s]", h->version, path);
  }

  if(strncmp(sig, "PATHS", 5) != 0) {
    if(!fatal) return -1;
    die("File is not valid paths file [path: %s]", path);
  }

  if(h->kmer_size % 2 == 0) {
    if(!fatal) return -1;
    die("kmer size is not an odd number [kmer_size: %u; binary: %s]\n",
        h->kmer_size, path);
  }

  if(h->kmer_size < 3) {
    if(!fatal) return -1;
    die("kmer size is less than three [kmer_size: %u; binary: %s]\n",
        h->kmer_size, path);
  }

  if(h->num_of_cols == 0) {
    if(!fatal) return -1;
    die("number of colours is zero [binary: %s]\n", path);
  }

  return bytes_read;
}

// We try loading the header into header passed
// Returns false if cannot read otherwise true
boolean paths_file_probe(const char *file_path, boolean *valid_paths_file,
                         PathFileHeader *pheader)
{
  FILE *fh = fopen(file_path, "r");
  if(fh == NULL) return false;
  int hret = paths_file_read_header(fh, pheader, false, file_path);
  fclose(fh);
  *valid_paths_file = (hret > 0);
  return true;
}

// Check a header and graph are compatible
void paths_graph_compatibility_check(const PathFileHeader *pheader,
                                     const dBGraph *db_graph)
{
  if(db_graph->kmer_size != pheader->kmer_size)
    die("Kmer sizes do not match between graph and path file");
  // if(db_graph->num_of_cols != pheader->num_of_cols)
  //   die("Number of colours does not match between graph and path file");

  if(pheader->num_path_bytes > db_graph->pdata.size) {
    die("Not enough memory allocated to store paths [mem: %zu]",
        (size_t)pheader->num_path_bytes);
  }

  if(db_graph->ht.unique_kmers > 0 &&
     db_graph->ht.unique_kmers < pheader->num_kmers_with_paths)
  {
    warn("Graph has fewer kmers than paths file");
  }

  // Check sample names match
  // uint32_t i;
  // for(i = 0; i < pheader->num_of_cols; i++)
  // {
  //   char *gname = db_graph->ginfo[i].sample_name.buff;
  //   char *pname = pheader->sample_names[i].buff;

  //   if(strcmp(pname, "noname") != 0 && strcmp(gname, pname) != 0)
  //     die("Graph/path sample names do not match [%u] '%s' vs '%s'", i, gname, pname);
  // }
}

// If tmppaths != NULL, do merge
// if insert is true, insert missing kmers into the graph
//  (this is a useful feature for pview)
void paths_format_merge(const char *path, PathFileHeader *pheader,
                        dBGraph *db_graph, PathStore *paths,
                        PathStore *tmppaths, boolean insert_missing_kmers)
{
  FILE *fh = fopen(path, "r");
  if(fh == NULL) die("Unable to open paths file: %s\n", path);
  setvbuf(fh, NULL, _IOFBF, CTP_BUF_SIZE);

  paths_file_read_header(fh, pheader, true, path);
  paths_graph_compatibility_check(pheader, db_graph);

  // Print some output
  char kmers_str[100], paths_str[100], mem_str[100];
  ulong_to_str(pheader->num_kmers_with_paths, kmers_str);
  ulong_to_str(pheader->num_of_paths, paths_str);
  bytes_to_str(pheader->num_path_bytes, 1, mem_str);

  status("Loading paths: %s paths, %s path-bytes, %s kmers\n",
         paths_str, mem_str, kmers_str);

  uint64_t i;
  BinaryKmer bkmer;
  hkey_t node;
  boolean found;

  // Load paths
  if(tmppaths != NULL)
  {
    if(pheader->num_path_bytes > tmppaths->size)
      die("Not enough memory for loading paths");
    safe_fread(fh, tmppaths->store, pheader->num_path_bytes, "tmppaths->store", path);
  }
  else
  {
    safe_fread(fh, paths->store, pheader->num_path_bytes, "paths->store", path);
    paths->next = paths->store + pheader->num_path_bytes;
    paths->num_of_paths = pheader->num_of_paths;
    paths->num_kmers_with_paths = pheader->num_kmers_with_paths;
  }

  // Load kmer pointers to paths
  PathIndex index;
  memset(bkmer.b, 0, sizeof(BinaryKmer));

  for(i = 0; i < pheader->num_kmers_with_paths; i++)
  {
    safe_fread(fh, &bkmer, sizeof(BinaryKmer), "bkmer", path);

    if(insert_missing_kmers) {
      node = hash_table_find_or_insert(&db_graph->ht, bkmer, &found);
    }
    else if((node = hash_table_find(&db_graph->ht, bkmer)) == HASH_NOT_FOUND) {
      die("Node missing: %zu [path: %s]", (size_t)node, path);
    }

    safe_fread(fh, &index, sizeof(uint64_t), "kmer_index", path);
    if(index > pheader->num_path_bytes) {
      die("Path index out of bounds [%zu > %zu]",
          (size_t)index, (size_t)pheader->num_path_bytes);
    }

    if(tmppaths == NULL)
      db_node_paths(db_graph, node) = index;
    else
    {
      // Do merge
      PathIndex kindex = db_node_paths(db_graph, node);
      do
      {
        kindex = path_store_add2(paths, kindex, tmppaths->store + index);
        if(kindex != PATH_NULL) db_node_paths(db_graph, node) = kindex;
        memcpy(&index, tmppaths->store + index, sizeof(PathIndex));
      }
      while(index != PATH_NULL);
    }
  }

  // Test that this is the end of the file
  uint8_t end;
  if(fread(&end, 1, 1, fh) != 0)
    warn("End of file not reached when loading! [path: %s]", path);

  fclose(fh);
}

void paths_format_read(const char *path, PathFileHeader *pheader,
                       dBGraph *db_graph, PathStore *paths,
                       boolean insert_missing_kmers)
{
  paths_format_merge(path, pheader, db_graph, paths, NULL, insert_missing_kmers);
}

//
// Write
//

// returns number of bytes written
size_t paths_format_write_header_core(const PathFileHeader *header, FILE *fout)
{
  fwrite("PATHS", 1, 5, fout);
  fwrite(&header->version, sizeof(uint32_t), 1, fout);
  fwrite(&header->kmer_size, sizeof(uint32_t), 1, fout);
  fwrite(&header->num_of_cols, sizeof(uint32_t), 1, fout);
  fwrite(&header->num_of_paths, sizeof(uint64_t), 1, fout);
  fwrite(&header->num_path_bytes, sizeof(uint64_t), 1, fout);
  fwrite(&header->num_kmers_with_paths, sizeof(uint64_t), 1, fout);
  return 5 + sizeof(uint32_t)*3 + sizeof(uint64_t)*3;
}

// returns number of bytes written
size_t paths_format_write_header(const PathFileHeader *header, FILE *fout)
{
  paths_format_write_header_core(header, fout);

  size_t i, bytes = 0;
  uint32_t len;
  const StrBuf *buf;

  for(i = 0; i < header->num_of_cols; i++)
  {
    buf = header->sample_names + i;
    len = buf->len;
    fwrite(&len, sizeof(uint32_t), 1, fout);
    fwrite(buf->buff, sizeof(uint8_t), len, fout);
    bytes += sizeof(uint32_t) + len;
  }

  return bytes;
}

static inline void write_optimised_paths(hkey_t node, PathIndex *pidx,
                                         dBGraph *db_graph, FILE *fout)
{
  PathStore *paths = &db_graph->pdata;
  PathIndex curridx, nextidx, newidx;
  PathLen len;
  Orientation orient;
  size_t mem;

  if((curridx = db_node_paths(db_graph, node)) != PATH_NULL)
  {
    db_node_paths(db_graph, node) = *pidx;

    do
    {
      nextidx = path_store_prev(paths, curridx);
      path_store_len_orient(paths, curridx, &len, &orient);
      mem = path_mem(paths->col_bitset_bytes,len);
      *pidx += mem;
      newidx = (nextidx == PATH_NULL ? PATH_NULL : *pidx);
      fwrite(&newidx, sizeof(PathIndex), 1, fout);
      fwrite(paths->store+curridx+sizeof(PathIndex), mem-sizeof(PathIndex), 1, fout);
      curridx = nextidx;
    }
    while(curridx != PATH_NULL);
  }
}

static inline void write_kmer_path_indices(hkey_t node, const dBGraph *db_graph,
                                           FILE *fout)
{
  if(db_node_paths(db_graph, node) != PATH_NULL)
  {
    BinaryKmer bkmer = db_node_bkmer(db_graph, node);
    PathIndex index = db_node_paths(db_graph, node);
    fwrite(&bkmer, sizeof(BinaryKmer), 1, fout);
    fwrite(&index, sizeof(PathIndex), 1, fout);
  }
}

// Corrupts paths so they cannot be used elsewhere
void paths_format_write_optimised_paths(dBGraph *db_graph, FILE *fout)
{
  PathIndex poffset = 0;
  HASH_TRAVERSE(&db_graph->ht, write_optimised_paths, &poffset, db_graph, fout);
  HASH_TRAVERSE(&db_graph->ht, write_kmer_path_indices, db_graph, fout);
}
