#include "global.h"
#include "gpath_reader.h"
#include "file_util.h"
#include "util.h"
#include "hash_mem.h"
#include "common_buffers.h"
#include "binary_seq.h"
#include "binary_kmer.h"
#include "gpath_set.h"
#include "gpath_store.h"
#include "gpath_subset.h"
#include "json_hdr.h"

/*
// File format:
<JSON_HEADER>
kmer [num] .. ignored
[FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. ignored
*/

#define load_check(x,msg,...) if(!(x)) { die("[LoadPathError] "msg, ##__VA_ARGS__); }

static long _json_demand_int(cJSON *root, const char *field, const char *path)
{
  cJSON *json = cJSON_GetObjectItem(root, field);
  if(json == NULL || json->type != cJSON_Number)
    die("No '%s' field in header: %s", field, path);
  return json->valueint;
}

size_t gpath_reader_get_kmer_size(const GPathReader *file)
{
  long val = _json_demand_int(file->json, "kmer_size", file->fltr.path.b);
  if(val < MIN_KMER_SIZE || val > MAX_KMER_SIZE || !(val & 1)) {
    die("kmer size is not an odd int between %i..%i: %li",
        MIN_KMER_SIZE, MAX_KMER_SIZE, val);
  }
  return val;
}

size_t gpath_reader_get_num_kmers(const GPathReader *file)
{
  // May be "num_kmers" or "kmers_with_paths"
  cJSON *json = cJSON_GetObjectItem(file->json, "num_kmers");
  if(json == NULL) json = cJSON_GetObjectItem(file->json, "num_kmers_with_paths");
  if(json == NULL || json->type != cJSON_Number) {
    die("No 'num_kmers' or 'num_kmers_with_paths' field in header: %s",
        file->fltr.path.b);
  }
  long val = json->valueint;
  // long val = _json_demand_int(file->json, "kmers_with_paths", file->fltr.path.b);
  if(val < 0) die("num_kmers is negative");
  return val;
}

size_t gpath_reader_get_num_paths(const GPathReader *file)
{
  long val = _json_demand_int(file->json, "num_paths", file->fltr.path.b);
  if(val < 0) die("num_paths is negative");
  return val;
}

size_t gpath_reader_get_path_bytes(const GPathReader *file)
{
  long val = _json_demand_int(file->json, "path_bytes", file->fltr.path.b);
  if(val < 0) die("path_bytes is negative");
  return val;
}

static size_t _gpath_reader_get_filencols(const GPathReader *file)
{
  long val = _json_demand_int(file->json, "ncols", file->fltr.path.b);
  if(val < 1 || val > 100000) die("Invalid number of colours: %li", val);
  return val;
}

// idx is 0..filencols
const char* gpath_reader_get_sample_name(const GPathReader *file, size_t idx)
{
  if(idx >= file->ncolours)
    die("Sample %zu > %zu: %s", idx, file->ncolours, file->fltr.path.b);
  cJSON *json = cJSON_GetObjectItem(file->colours_json[idx], "sample");
  if(json == NULL || json->type != cJSON_String)
    die("Sample %zu has no field 'sample': %s", idx, file->fltr.path.b);
  return json->valuestring;
}

void _parse_json_header(GPathReader *file)
{
  cJSON *sample, *colours = cJSON_GetObjectItem(file->json, "colours");
  if(colours == NULL || colours->type != cJSON_Array)
    die("No 'colours' array entry in header");

  file->ncolours = cJSON_GetArraySize(colours);
  file->colours_json = ctx_calloc(file->ncolours, sizeof(cJSON*));
  size_t i = 0;
  for(sample = colours->child; sample != NULL; sample = sample->next)
    file->colours_json[i++] = sample;

  if(file->ncolours == 0) die("No colours in JSON header");
}

// Open file, exit on error
// if successful creates a new GPathReader and returns 1
void gpath_reader_open2(GPathReader *file, char *path, const char *mode)
{
  FileFilter *fltr = &file->fltr;
  file_filter_open(fltr, path); // calls die() on error

  file->gz = futil_gzopen(fltr->path.b, mode);

  // Load JSON header into file->hdrstr
  StrBuf *hdrstr = &file->hdrstr;
  if(hdrstr->b == NULL) strbuf_alloc(hdrstr, 1024);
  json_hdr_read(NULL, file->gz, path, hdrstr);
  file->json = cJSON_Parse(hdrstr->b);
  if(file->json == NULL) die("Invalid JSON header: %s", path);

  // DEV: validate json header
  _parse_json_header(file);

  // The following functions call die() if kmer_size or ncols header fields
  // are missing
  size_t kmer_size = gpath_reader_get_kmer_size(file);
  size_t filencols = _gpath_reader_get_filencols(file);
  file_filter_set_cols(fltr, filencols);

  // Check we can handle the kmer size
  db_graph_check_kmer_size(kmer_size, file->fltr.path.b);
}

void gpath_reader_open(GPathReader *file, char *path)
{
  gpath_reader_open2(file, path, "r");
}

void gpath_reader_close(GPathReader *file)
{
  if(file->gz) gzclose(file->gz);
  file_filter_close(&file->fltr);
  cJSON_Delete(file->json);
  strbuf_dealloc(&file->hdrstr);
  ctx_free(file->colours_json);
  memset(file, 0, sizeof(GPathReader));
}


void gpath_reader_check(const GPathReader *file, size_t db_kmer_size,
                        size_t db_ncols)
{
  size_t kmer_size = gpath_reader_get_kmer_size(file);
  if(kmer_size != db_kmer_size) {
    die("Path file kmer size doesn't match [path: %s kmer: %zu, graph-kmer: %zu]",
        file->fltr.path.b, kmer_size, db_kmer_size);
  }
  size_t used_ncols = file_filter_into_ncols(&file->fltr);
  if(used_ncols > db_ncols) {
    die("Path file requires at least %zu colours [path: %s, curr colours: %zu]",
        used_ncols, file->fltr.input.b, db_ncols);
  }
}

// <KMER> <num> .. (ignored)
static void _gpath_reader_load_kmer_line(const char *path,
                                         StrBuf *line, dBGraph *db_graph,
                                         bool dont_add_kmers,
                                         size_t *npaths_ptr, hkey_t *hkey_ptr)
{
  const size_t kmer_size = db_graph->kmer_size;
  char *numpstr, *endpstr;
  size_t i, num_paths_exp = 0;

  for(i = 0; i < line->end && char_is_acgt(line->b[i]); i++) {}
  load_check(i == kmer_size, "Bad kmer line: %s\n'%s'\n", path, line->b);

  BinaryKmer bkmer, bkey;
  hkey_t hkey;
  bool found = false;

  bkmer = binary_kmer_from_str(line->b, kmer_size);
  // check bkmer is kmer key
  bkey = binary_kmer_get_key(bkmer, kmer_size);
  load_check(binary_kmers_are_equal(bkmer, bkey), "Bkmer not bkey: %s", path);

  if(dont_add_kmers) {
    hkey = hash_table_find(&db_graph->ht, bkey);
    load_check(hkey != HASH_NOT_FOUND, "BKmer not already loaded: %s", path);
  } else {
    hkey = hash_table_find_or_insert(&db_graph->ht, bkey, &found);
  }

  // Parse number of paths
  numpstr = line->b+kmer_size;
  load_check(numpstr[0] == ' ', "Bad kmer line: %s", path);
  numpstr++;
  num_paths_exp = strtoul(numpstr, &endpstr, 10);
  load_check(endpstr > numpstr && (*endpstr == '\0' || *endpstr == ' '),
              "Bad kmer line: %s\n%s\n", path, line->b);

  *npaths_ptr = num_paths_exp;
  *hkey_ptr = hkey;
}

// [FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. (ignored)
static void _gpath_reader_load_path_line(GPathReader *file, const char *path,
                                         StrBuf *line,
                                         uint8_t *nseenbuf, ByteBuffer *seqbuf,
                                         GPathSet *gpset)
{
  size_t i, num_kmers, num_juncs;
  Orientation orient;
  char *pstr, *endpstr;

  orient = (line->b[0] == 'F') ? FORWARD : REVERSE;
  load_check(line->b[1] == ' ', "Bad path line: %s", path);
  pstr = line->b+2;

  // [nkmers]
  num_kmers = strtoul(pstr, &endpstr, 10);
  load_check(endpstr > pstr && *endpstr == ' ', "Bad path line: %s", path);
  pstr = endpstr+1;

  // [njuncs]
  num_juncs = strtoul(pstr, &endpstr, 10);
  load_check(endpstr > pstr && *endpstr == ' ', "Bad path line: %s", path);
  pstr = endpstr+1;

  // [nseen,nseen,...]
  endpstr = strchr(pstr, ' ');
  load_check(endpstr != NULL && endpstr > pstr, "Bad path line: %s", path);

  for(i = 0; i < file->fltr.filencols; i++, pstr = endpstr+1)
  {
    size_t num_seen = strtoul(pstr, &endpstr, 10);
    nseenbuf[i] = MIN2(UINT8_MAX, num_seen);
    load_check((i+1 < file->fltr.filencols && *endpstr == ',') || *endpstr == ' ',
               "Bad path line: %s\n%s", path, line->b);
  }

  ctx_assert2(num_kmers < GPATH_MAX_KMERS, "%zu", (size_t)num_kmers);
  ctx_assert2(num_juncs < GPATH_MAX_JUNCS, "%zu", (size_t)num_juncs);

  // [seq]
  endpstr = pstr;
  do {
    load_check(char_is_acgt(*endpstr), "Bad path line: %s", path);
    endpstr++;
  } while(*endpstr && *endpstr != ' ');
  load_check((size_t)(endpstr - pstr) == num_juncs, "Bad path line: %s", path);
  byte_buf_ensure_capacity(seqbuf, (num_juncs+3)/4);
  binary_seq_from_str(pstr, num_juncs, seqbuf->data);

    // Add to GPathSet
  GPathNew newgpath = {.seq = seqbuf->data,
                       .colset = NULL, .nseen = NULL,
                       .orient = orient,
                       .klen = num_kmers,
                       .num_juncs = num_juncs};

  GPath *gpath = gpath_set_add_mt(gpset, newgpath);

  // Sort out colour re-ordering, and record nseen / colour bitset
  size_t from, to;
  uint8_t *nseen = gpath_set_get_nseen(gpset, gpath);
  uint8_t *colset = gpath_get_colset(gpath, gpset->ncols);

  memset(nseen, 0, gpset->ncols);

  for(i = 0; i < file->fltr.ncols; i++) {
    from = file_filter_fromcol(&file->fltr, i);
    to = file_filter_intocol(&file->fltr, i);
    nseen[to] = MIN2((size_t)UINT8_MAX, (size_t)nseen[to] + nseenbuf[from]);
    bitset_or(colset, to, nseen[to] > 0);
  }
}

static void _load_paths_from_set(dBGraph *db_graph, GPathSet *gpset,
                                 GPathSubset *subset0, GPathSubset *subset1,
                                 hkey_t hkey, size_t num_paths_exp,
                                 const char *file_path)
{
  GPathStore *gpstore = &db_graph->gpstore;
  GPathHash *gphash = &db_graph->gphash;

  size_t i;

  load_check(gpset->entries.len == num_paths_exp,
             "Too many/few paths: %s (exp %zu vs act %zu)",
             file_path, num_paths_exp, gpset->entries.len);

  gpath_subset_init(subset0, &gpstore->gpset);
  gpath_subset_init(subset1, gpset);

  GPath *kmer_paths = gpath_store_fetch(gpstore, hkey);
  gpath_subset_load_llist(subset0, kmer_paths);
  gpath_subset_load_set(subset1);
  gpath_subset_rmdup(subset1);

  // Merge entries from subset1 into subset0
  gpath_subset_merge(subset0, subset1, false);

  bool found;

  // Copy remaining entries across
  for(i = 0; i < subset1->list.len; i++)
  {
    GPathNew newgp = gpath_set_get(gpset, subset1->list.data[i]);

    if(db_graph_has_path_hash(db_graph))
      gpath_hash_find_or_insert_mt(gphash, hkey, newgp, &found);
    else
      gpath_store_add_mt(gpstore, hkey, newgp);
  }

  gpath_set_reset(gpset);
}

void gpath_reader_load(GPathReader *file, bool dont_add_kmers, dBGraph *db_graph)
{
  gzFile gzin = file->gz;
  const char *path = file->fltr.path.b;

  status("Loading path file: %s", path);

  // Read kmer lines
  // <KMER> <num>
  // [FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. ignored
  StrBuf line;
  strbuf_alloc(&line, 2048);

  uint8_t *nseenbuf = ctx_calloc(file->fltr.filencols, sizeof(uint8_t));

  ByteBuffer seqbuf;
  byte_buf_alloc(&seqbuf, 64);

  // Load paths into this temporary set for each kmer
  GPathSet gpset;
  gpath_set_alloc(&gpset, db_graph->num_of_cols, ONE_MEGABYTE, true, true);

  GPathSubset subset0, subset1;
  gpath_subset_alloc(&subset0);
  gpath_subset_alloc(&subset1);

  hkey_t hkey = HASH_NOT_FOUND;
  size_t num_paths_exp = 0, num_kmers = 0;

  size_t num_kmers_exp = gpath_reader_get_num_kmers(file);

  while(strbuf_reset_gzreadline(&line, gzin) > 0) {
    strbuf_chomp(&line);
    if(line.end > 0 && line.b[0] != '#') {
      if(line.b[0] == 'F' || line.b[0] == 'R') {
        // Assume path line
        // status("Path line: %s", line.b);
        load_check(hkey != HASH_NOT_FOUND, "Path before kmer: %s", path);

        _gpath_reader_load_path_line(file, path, &line, nseenbuf, &seqbuf,
                                     &gpset);

        load_check(gpset.entries.len <= num_paths_exp, "Too many paths: %s", path);
      }
      else {
        // Assume kmer line
        // status("Kmer line: %s", line.b);

        if(num_kmers > 0) {
          _load_paths_from_set(db_graph, &gpset, &subset0, &subset1,
                               hkey, num_paths_exp, path);
        }

        _gpath_reader_load_kmer_line(path, &line, db_graph, dont_add_kmers,
                                     &num_paths_exp, &hkey);

        num_kmers++;
      }
    }
  }

  if(num_kmers > 0) {
    _load_paths_from_set(db_graph, &gpset, &subset0, &subset1,
                         hkey, num_paths_exp, path);
  }

  load_check(num_kmers == (size_t)num_kmers_exp,
             "num_kmers don't match (exp %zu vs %zu)",
             num_kmers, (size_t)num_kmers_exp);

  gpath_subset_dealloc(&subset0);
  gpath_subset_dealloc(&subset1);

  gpath_set_dealloc(&gpset);
  strbuf_dealloc(&line);
  byte_buf_dealloc(&seqbuf);
  ctx_free(nseenbuf);
}

void gpath_reader_load_sample_names(const GPathReader *file, dBGraph *db_graph)
{
  const FileFilter *fltr = &file->fltr;

  size_t i, intocol, fromcol;
  StrBuf *gname;
  const char *pname;

  for(i = 0; i < fltr->ncols; i++)
  {
    intocol = file_filter_intocol(fltr, i);
    fromcol = file_filter_fromcol(fltr, i);
    gname = &db_graph->ginfo[intocol].sample_name;
    pname = gpath_reader_get_sample_name(file, fromcol);

    if(strcmp(gname->b, "undefined") == 0) {
      strbuf_set(gname, pname);
    } else if(strcmp(gname->b, pname) != 0) {
      die("Graph/path sample names do not match [%zu->%zu] '%s' vs '%s'",
          fromcol, intocol, gname->b, pname);
    }
  }
}

//
// Memory Calculations
//

// Get max mem required to load ctp files
// Assume equal three-way split between: GPath, seq+colset, hashtable
void gpath_reader_max_mem_req(GPathReader *files, size_t nfiles,
                              size_t ncols, size_t graph_capacity,
                              bool store_nseen_klen,
                              bool split_lists, bool use_hash,
                              size_t *min_mem_ptr, size_t *max_mem_ptr)
{
  size_t i, path_bytes, hash_bytes;
  size_t max_npaths = 0, sum_npaths = 0, max_pbytes = 0, sum_pbytes = 0;

  path_bytes = sizeof(GPath) +
               (store_nseen_klen ? sizeof(uint8_t)*ncols + sizeof(uint32_t) : 0);
  hash_bytes = (use_hash ? sizeof(GPEntry)/IDEAL_OCCUPANCY : 0);

  for(i = 0; i < nfiles; i++) {
    size_t npaths = gpath_reader_get_num_paths(&files[i]);
    size_t pbytes = gpath_reader_get_path_bytes(&files[i]) + npaths*(ncols+7)/8;
    max_npaths = MAX2(max_npaths, npaths);
    max_pbytes = MAX2(max_pbytes, pbytes);
    sum_npaths += npaths;
    sum_pbytes += pbytes;
  }

  size_t list_mem = graph_capacity*sizeof(GPath*) * (split_lists ? 2 : 1);

  // Memory is split three ways between path entries, sequence+colsets, hashtable
  size_t mult = use_hash ? 3 : 2;
  *min_mem_ptr = MAX3(max_npaths*path_bytes, max_pbytes, max_npaths*hash_bytes)*mult
                 + list_mem;
  *max_mem_ptr = MAX3(sum_npaths*path_bytes, sum_pbytes, sum_npaths*hash_bytes)*mult
                 + list_mem;
}

size_t gpath_reader_sum_mem(GPathReader *files, size_t nfiles,
                            size_t ncols, bool count_nseen, bool use_gphash,
                            size_t *max_file_mem_ptr)
{
  size_t i, npaths, path_bytes, file_mem;
  size_t path_sum_mem = 0, max_file_mem = 0;

  for(i = 0; i < nfiles; i++)
  {
    npaths = gpath_reader_get_num_paths(&files[i]);
    path_bytes = gpath_reader_get_path_bytes(&files[i]);
    file_mem = npaths * sizeof(GPath) + // GPath
               path_bytes + // Sequence
               npaths * (ncols+7)/8 + // Colset
               npaths * (count_nseen ? ncols*sizeof(uint8_t) + sizeof(uint32_t) : 0) +
               (npaths/IDEAL_OCCUPANCY) * (use_gphash ? sizeof(GPEntry) : 0);

    path_sum_mem += file_mem;
    max_file_mem = MAX2(max_file_mem, file_mem);
  }

  if(max_file_mem_ptr) *max_file_mem_ptr = max_file_mem;

  return path_sum_mem;
}

size_t gpath_reader_mem_req(GPathReader *files, size_t nfiles,
                            size_t ncols, size_t max_mem,
                            bool count_nseen)
{
  size_t max_file_mem = 0, path_sum_mem;
  path_sum_mem = gpath_reader_sum_mem(files, nfiles, ncols, count_nseen, false,
                                      &max_file_mem);

  if(max_file_mem > max_mem) {
    char memstr[50];
    bytes_to_str(max_file_mem, 1, memstr);
    die("Require at least %s memory for paths", memstr);
  }

  return MIN2(max_mem, path_sum_mem);
}

// Create a path store that does not tracks path counts
void gpath_reader_alloc_gpstore(GPathReader *files, size_t nfiles,
                                size_t mem, bool count_nseen,
                                dBGraph *db_graph)
{
  if(nfiles == 0) return;

  size_t sum_mem = gpath_reader_sum_mem(files, nfiles, db_graph->num_of_cols,
                                        count_nseen, false, NULL);

  size_t i, npaths, max_npaths = 0, sum_npaths = 0;

  for(i = 0; i < nfiles; i++) {
    npaths = gpath_reader_get_num_paths(&files[i]);
    max_npaths = MAX2(max_npaths, npaths);
    sum_npaths += npaths;
  }

  npaths = sum_mem <= mem ? sum_npaths : 0;
  status("[GPathReader] need %zu paths", npaths);

  gpath_store_alloc(&db_graph->gpstore,
                    db_graph->num_of_cols,
                    db_graph->ht.capacity,
                    npaths, mem,
                    count_nseen, false);
}
