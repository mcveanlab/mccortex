#include "global.h"
#include "path_format.h"
#include "db_graph.h"
#include "db_node.h"
#include "hash_table.h"
#include "binary_paths.h"
#include "util.h"
#include "file_util.h"

// Format:
// "PATHS"<uint32_t:version><uint32_t:kmersize><uint32_t:num_cols>
// <uint64_t:num_of_paths><uint64_t:num_path_bytes><uint64_t:num_path_kmers>
// <uint8_t:path_data>
// <binarykmer><uint64_t:path_index_fw><uint64_t:path_index_rv>

typedef struct
{
  uint32_t version, kmer_size, num_cols;
  uint64_t num_of_paths, num_path_bytes, num_path_kmers;
} PathFileHeader;

static inline void write_kmer_path_indices(hkey_t node, const dBGraph *db_graph,
                                           FILE *fout)
{
  if(db_node_paths(db_graph, node) != PATH_NULL)
  {
    fwrite(db_node_bkmer(db_graph, node), sizeof(BinaryKmer), 1, fout);
    uint64_t index = db_node_paths(db_graph, node);
    fwrite(&index, sizeof(uint64_t), 1, fout);
  }
}

void paths_format_write(const dBGraph *db_graph, const PathStore *paths,
                        const char *path)
{
  FILE *fout;

  if((fout = fopen(path, "w")) == NULL) {
    die("Unable to open binary paths file to write: %s\n", path);
  }

  fwrite("PATHS", 1, 5, fout);

  uint32_t version = CTX_PATH_FILEFORMAT;
  uint32_t kmer_size = db_graph->kmer_size;
  uint32_t num_cols = db_graph->num_of_cols_used;
  fwrite(&version, sizeof(uint32_t), 1, fout);
  fwrite(&kmer_size, sizeof(uint32_t), 1, fout);
  fwrite(&num_cols, sizeof(uint32_t), 1, fout);

  uint64_t num_of_paths = paths->num_of_paths;
  uint64_t num_path_bytes = paths->next - paths->store;
  uint64_t num_kmers_with_paths = paths->num_kmers_with_paths;

  char kmers_str[100], paths_str[100], mem_str[100];
  ulong_to_str(num_kmers_with_paths, kmers_str);
  ulong_to_str(num_of_paths, paths_str);
  bytes_to_str(num_path_bytes, 1, mem_str);

  message("  Saving paths: %s paths, %s path-bytes, %s kmers\n",
          paths_str, mem_str, kmers_str);

  fwrite(&num_of_paths, sizeof(uint64_t), 1, fout);
  fwrite(&num_path_bytes, sizeof(uint64_t), 1, fout);
  fwrite(&num_kmers_with_paths, sizeof(uint64_t), 1, fout);
  fwrite(paths->store, sizeof(uint8_t), num_path_bytes, fout);

  HASH_TRAVERSE(&db_graph->ht, write_kmer_path_indices, db_graph, fout);

  fclose(fout);
}

static void paths_format_load_header(PathFileHeader *h,
                                     FILE *fh, const char *path,
                                     PathStore *paths, dBGraph *db_graph)
{
  char sig[6] = {0};

  safe_fread(fh, sig, 5, "PATHS", path);
  safe_fread(fh, &h->version, sizeof(uint32_t), "version", path);
  safe_fread(fh, &h->kmer_size, sizeof(uint32_t), "kmer_size", path);
  safe_fread(fh, &h->num_cols, sizeof(uint32_t), "num_cols", path);
  safe_fread(fh, &h->num_of_paths, sizeof(uint64_t), "num_of_paths", path);
  safe_fread(fh, &h->num_path_bytes, sizeof(uint64_t), "num_path_bytes", path);
  safe_fread(fh, &h->num_path_kmers, sizeof(uint64_t), "num_path_kmers", path);

  if(strncmp(sig, "PATHS", 5) != 0) {
    die("File is not valid paths file [path: %s]", path);
  }
  if(h->version != CTX_PATH_FILEFORMAT) {
    die("file version not supported [version: %u; path: %s]", h->version, path);
  }
  if(h->kmer_size != db_graph->kmer_size) {
    die("kmer_size values don't match [%u vs %u; path: %s]",
        h->kmer_size, db_graph->kmer_size, path);
  }
  if(h->num_cols != paths->num_of_cols) {
    die("numbers of colours don't match [%u vs %u; path: %s]",
        h->num_cols, paths->num_of_cols, path);
  }
  if(h->num_path_bytes > paths->size) {
    die("Not enough memory allocated to store paths from file: %s [mem: %zu]",
        path, (size_t)h->num_path_bytes);
  }
}

// If tmppaths != NULL, do merge
// if insert is true, insert missing kmers into the graph
//  (this is a useful feature for pview)
void paths_format_read(dBGraph *db_graph, PathStore *paths, PathStore *tmppaths,
                       boolean insert_missing_kmers, const char *path)
{
  FILE *fh;

  if((fh = fopen(path, "r")) == NULL) {
    die("Unable to open binary paths file to read: %s\n", path);
  }

  PathFileHeader header;
  paths_format_load_header(&header, fh, path, paths, db_graph);

  if(header.kmer_size != db_graph->kmer_size) {
    die("kmer_size values don't match [%u vs %u; path: %s]",
        header.kmer_size, db_graph->kmer_size, path);
  }

  char kmers_str[100], paths_str[100], mem_str[100];
  ulong_to_str(header.num_path_kmers, kmers_str);
  ulong_to_str(header.num_of_paths, paths_str);
  bytes_to_str(header.num_path_bytes, 1, mem_str);

  message(" Loading paths: %s paths, %s path-bytes, %s kmers\n",
          paths_str, mem_str, kmers_str);

  if(header.num_path_bytes > tmppaths->size)
    die("Not enough memory for loading paths");

  if(tmppaths != NULL)
    safe_fread(fh, tmppaths->store, header.num_path_bytes, "tmppaths->store", path);
  else
  {
    safe_fread(fh, paths->store, header.num_path_bytes, "paths->store", path);
    paths->next = paths->store + header.num_path_bytes;
    paths->num_of_paths = header.num_of_paths;
    paths->num_kmers_with_paths = header.num_path_kmers;
  }

  BinaryKmer bkmer;
  uint64_t i;
  PathIndex index;
  hkey_t node;
  boolean found;
  memset(bkmer, 0, sizeof(BinaryKmer));

  for(i = 0; i < header.num_path_kmers; i++)
  {
    safe_fread(fh, &bkmer, sizeof(BinaryKmer), "bkmer", path);

    if(insert_missing_kmers) {
      node = hash_table_find_or_insert(&db_graph->ht, bkmer, &found);
    }
    else if((node = hash_table_find(&db_graph->ht, bkmer)) == HASH_NOT_FOUND) {
      die("Node missing: %zu [path: %s]", (size_t)node, path);
    }

    safe_fread(fh, &index, sizeof(uint64_t), "kmer_index", path);
    if(index > header.num_path_bytes) {
      die("Path index out of bounds [%zu > %zu]",
          (size_t)index, (size_t)header.num_path_bytes);
    }

    if(tmppaths == NULL)
      db_node_paths(db_graph, node) = index;
    else
    {
      // Do merge
      PathIndex kindex = db_node_paths(db_graph, node);
      do
      {
        kindex = binary_paths_add2(paths, kindex, tmppaths->store + index);
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

// Returns false if cannot read otherwise true
boolean paths_format_probe(const char *path, boolean *valid_paths_file,
                           uint32_t *kmer_size_ptr, uint32_t *num_of_cols_ptr,
                           uint64_t *num_of_paths_ptr, uint64_t *num_path_bytes_ptr,
                           uint64_t *num_path_kmers_ptr)
{
  FILE *fh;

  if((fh = fopen(path, "r")) == NULL) return false;

  char header[6] = {0};
  uint32_t version, kmer_size, num_cols;
  uint64_t num_of_paths, num_path_bytes, num_path_kmers;

  *valid_paths_file = false;

  if(fread(header, 1, 5, fh) == 5 &&
     fread(&version, sizeof(uint32_t), 1, fh) == 1 &&
     fread(&kmer_size, sizeof(uint32_t), 1, fh) == 1 &&
     fread(&num_cols, sizeof(uint32_t), 1, fh) == 1 &&
     fread(&num_of_paths, sizeof(uint64_t), 1, fh) == 1 &&
     fread(&num_path_bytes, sizeof(uint64_t), 1, fh) == 1 &&
     fread(&num_path_kmers, sizeof(uint64_t), 1, fh) == 1 &&
     strncmp(header, "PATHS", 5) == 0 && version == CTX_PATH_FILEFORMAT)
  {
    *valid_paths_file = true;
    *kmer_size_ptr = kmer_size;
    *num_of_cols_ptr = num_cols;
    *num_of_paths_ptr = num_of_paths;
    *num_path_bytes_ptr = num_path_bytes;
    *num_path_kmers_ptr = num_path_kmers;
  }

  // printf("header:%s, version:%u, kmer_size:%u, num_cols:%u\n",
  //        header, version, kmer_size, num_cols);

  fclose(fh);
  return true;
}

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
