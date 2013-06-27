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
// <uint64_t:num_paths><uint64_t:num_path_bytes><uint64_t:num_path_kmers>
// <uint8_t:path_data>
// <binarykmer><uint64_t:path_index_fw><uint64_t:path_index_rv>

static inline void count_kmers_with_paths(hkey_t node, const dBGraph *db_graph,
                                          uint64_t *count)
{
  if(db_node_paths(db_graph, node, forward) != PATH_NULL ||
     db_node_paths(db_graph, node, reverse) != PATH_NULL)
  {
    (*count)++;
  }
}

static inline void write_kmer_path_indices(hkey_t node, const dBGraph *db_graph,
                                           FILE *fout)
{
  if(db_node_paths(db_graph, node, forward) != PATH_NULL ||
     db_node_paths(db_graph, node, reverse) != PATH_NULL)
  {
    fwrite(db_node_bkmer(db_graph, node), sizeof(BinaryKmer), 1, fout);
    fwrite(&db_node_paths(db_graph, node, forward), sizeof(uint64_t), 1, fout);
    fwrite(&db_node_paths(db_graph, node, reverse), sizeof(uint64_t), 1, fout);
  }
}

void paths_format_write(const dBGraph *db_graph, const binary_paths_t *paths,
                        const char *path)
{
  FILE *fout;

  if((fout = fopen(path, "w")) == NULL) {
    die("Unable to open binary paths file to write: %s\n", path);
  }

  fwrite("PATHS", 1, 5, fout);

  uint32_t version = CTX_PATH_FILEFORMAT;
  uint32_t kmer_size = db_graph->kmer_size;
  uint32_t num_cols = db_graph->ginfo.num_of_colours_loaded;
  fwrite(&version, sizeof(uint32_t), 1, fout);
  fwrite(&kmer_size, sizeof(uint32_t), 1, fout);
  fwrite(&num_cols, sizeof(uint32_t), 1, fout);

  uint64_t num_paths = paths->num_paths;
  uint64_t num_path_bytes = paths->next - paths->store;
  uint64_t num_path_kmers = 0;
  HASH_TRAVERSE(&db_graph->ht, count_kmers_with_paths, db_graph, &num_path_kmers);

  char kmers_str[100], paths_str[100], mem_str[100];
  ulong_to_str(num_path_kmers, kmers_str);
  ulong_to_str(num_paths, paths_str);
  bytes_to_str(num_path_bytes, 1, mem_str);

  message("  Saving %s kmers holding %s paths (%s)\n",
          kmers_str, paths_str, mem_str);

  fwrite(&num_paths, sizeof(uint64_t), 1, fout);
  fwrite(&num_path_bytes, sizeof(uint64_t), 1, fout);
  fwrite(&num_path_kmers, sizeof(uint64_t), 1, fout);
  fwrite(paths->store, sizeof(uint8_t), num_path_bytes, fout);

  HASH_TRAVERSE(&db_graph->ht, write_kmer_path_indices, db_graph, fout);

  fclose(fout);
}

// if insert is true, insert missing kmers into the graph
void paths_format_read(dBGraph *db_graph, binary_paths_t *paths,
                       boolean insert, const char *path)
{
  FILE *fh;

  if((fh = fopen(path, "r")) == NULL) {
    die("Unable to open binary paths file to read: %s\n", path);
  }

  char header[6] = {0};
  uint32_t version, kmer_size = db_graph->kmer_size, num_cols;
  uint64_t num_paths, num_path_bytes, num_path_kmers;

  safe_fread(fh, header, 5, "PATHS", path);
  safe_fread(fh, &version, sizeof(uint32_t), "version", path);
  safe_fread(fh, &kmer_size, sizeof(uint32_t), "kmer_size", path);
  safe_fread(fh, &num_cols, sizeof(uint32_t), "num_cols", path);
  safe_fread(fh, &num_paths, sizeof(uint64_t), "num_paths", path);
  safe_fread(fh, &num_path_bytes, sizeof(uint64_t), "num_path_bytes", path);
  safe_fread(fh, &num_path_kmers, sizeof(uint64_t), "num_path_kmers", path);

  if(strncmp(header, "PATHS", 5) != 0) {
    die("File is not valid paths file [path: %s]", path);
  }
  if(version != CTX_PATH_FILEFORMAT) {
    die("file version not supported [version: %u; path: %s]", version, path);
  }
  if(kmer_size != db_graph->kmer_size) {
    die("kmer_size values don't match [%u vs %u; path: %s]",
        kmer_size, db_graph->kmer_size, path);
  }
  if(num_cols != db_graph->ginfo.num_of_colours_loaded) {
    die("numbers of colours don't match [%u vs %u; path: %s]",
        num_cols, db_graph->ginfo.num_of_colours_loaded, path);
  }
  if(num_path_bytes > paths->size) {
    die("Not enough memory allocated to store paths from file: %s [mem: %zu]",
        path, (size_t)num_path_bytes);
  }

  message(" Loaded paths: %zu paths, %zu path-bytes, %zu kmers\n",
          (size_t)num_paths, (size_t)num_path_bytes, (size_t)num_path_kmers);

  safe_fread(fh, paths->store, num_path_bytes, "paths->store", path);

  BinaryKmer bkmer;
  uint64_t i, index;
  hkey_t node;
  boolean found;
  memset(bkmer, 0, sizeof(BinaryKmer));

  for(i = 0; i < num_path_kmers; i++)
  {
    safe_fread(fh, &bkmer, sizeof(BinaryKmer), "bkmer", path);

    if(insert) {
      node = hash_table_find_or_insert(&db_graph->ht, bkmer, &found);
    }
    else if((node = hash_table_find(&db_graph->ht, bkmer)) == HASH_NOT_FOUND) {
      die("Node missing: %zu [path: %s]", (size_t)node, path);
    }

    safe_fread(fh, &index, sizeof(uint64_t), "fw_index", path);
    db_node_paths(db_graph, node, forward) = index;

    safe_fread(fh, &index, sizeof(uint64_t), "rv_index", path);
    db_node_paths(db_graph, node, reverse) = index;
  }

  paths->next = paths->store + num_path_bytes;
  paths->num_paths = num_paths;

  // Test that this is the end of the file
  uint8_t end;
  if(fread(&end, 1, 1, fh) != 0)
    warn("End of file not reached when loading! [path: %s]", path);

  fclose(fh);
}

// Returns false if cannot read otherwise true
boolean paths_format_probe(const char *path, boolean *valid_paths_file,
                           uint32_t *kmer_size_ptr, uint32_t *num_of_cols_ptr,
                           uint64_t *num_paths_ptr, uint64_t *num_path_bytes_ptr,
                           uint64_t *num_path_kmers_ptr)
{
  FILE *fh;

  if((fh = fopen(path, "r")) == NULL) return false;

  char header[6] = {0};
  uint32_t version, kmer_size, num_cols;
  uint64_t num_paths, num_path_bytes, num_path_kmers;

  *valid_paths_file = false;

  if(fread(header, 1, 5, fh) == 5 &&
     fread(&version, sizeof(uint32_t), 1, fh) == 1 &&
     fread(&kmer_size, sizeof(uint32_t), 1, fh) == 1 &&
     fread(&num_cols, sizeof(uint32_t), 1, fh) == 1 &&
     fread(&num_paths, sizeof(uint64_t), 1, fh) == 1 &&
     fread(&num_path_bytes, sizeof(uint64_t), 1, fh) == 1 &&
     fread(&num_path_kmers, sizeof(uint64_t), 1, fh) == 1 &&
     strncmp(header, "PATHS", 5) == 0 && version == CTX_PATH_FILEFORMAT)
  {
    *valid_paths_file = true;
    *kmer_size_ptr = kmer_size;
    *num_of_cols_ptr = num_cols;
    *num_paths_ptr = num_paths;
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
  if(strcasecmp(path+len-4,".ctx") != 0) len += 4;
  memcpy(out+len-4, ".ctp", 4);
  out[len] = '\0';
}
