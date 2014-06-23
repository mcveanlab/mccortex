#include "global.h"
#include "gpath_reader.h"
#include "util.h"
#include "hash_mem.h"
#include "common_buffers.h"
#include "binary_seq.h"
#include "binary_kmer.h"
#include "gpath_set.h"
#include "gpath_store.h"
#include "gpath_subset.h"

#define MAX_JSON_HDR_BYTES (1<<20) /* 1M max json header */

/*
// File format:
<JSON_HEADER>
kmer [num] .. ignored
[FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. ignored
*/

#define load_assert(x,msg,...) if(!(x)) { die("[LoadPathError] "msg, ##__VA_ARGS__); }

static int _json_demand_int(cJSON *root, const char *field, const char *path)
{
  cJSON *json = cJSON_GetObjectItem(root, field);
  if(json == NULL || json->type != cJSON_Number)
    die("No '%s' field in header: %s", field, path);
  return json->valueint;
}

size_t gpath_reader_get_kmer_size(const GPathReader *file)
{
  int val = _json_demand_int(file->json, "kmer_size", file->fltr.file_path.buff);
  if(val < MIN_KMER_SIZE || val > MAX_KMER_SIZE || !(val & 1)) {
    die("kmer size is not an odd int between %i..%i: %i",
        MIN_KMER_SIZE, MAX_KMER_SIZE, val);
  }
  return val;
}

size_t gpath_reader_get_num_kmers(const GPathReader *file)
{
  int val = _json_demand_int(file->json, "num_kmers", file->fltr.file_path.buff);
  if(val < 0) die("num_kmers is negative");
  return val;
}

size_t gpath_reader_get_num_paths(const GPathReader *file)
{
  int val = _json_demand_int(file->json, "num_paths", file->fltr.file_path.buff);
  if(val < 0) die("num_paths is negative");
  return val;
}

size_t gpath_reader_get_path_bytes(const GPathReader *file)
{
  int val = _json_demand_int(file->json, "path_bytes", file->fltr.file_path.buff);
  if(val < 0) die("path_bytes is negative");
  return val;
}

static size_t _gpath_reader_get_filencols(const GPathReader *file)
{
  int val = _json_demand_int(file->json, "ncols", file->fltr.file_path.buff);
  if(val < 1 || val > 100000) die("Invalid number of colours: %i", val);
  return val;
}

// idx is 0..filencols
const char* gpath_reader_get_sample_name(const GPathReader *file, size_t idx)
{
  if(idx >= file->ncolours)
    die("Sample %zu > %zu: %s", idx, file->ncolours, file->fltr.file_path.buff);
  cJSON *json = cJSON_GetObjectItem(file->colours_json[idx], "sample");
  if(json == NULL || json->type != cJSON_String)
    die("Sample %zu has no field 'sample': %s", idx, file->fltr.file_path.buff);
  return json->valuestring;
}

static void _gpath_reader_load_json_hdr(gzFile gzin, const char *path,
                                        StrBuf *hdrstr)
{
  size_t nread;
  strbuf_reset(hdrstr);

  load_assert((nread = strbuf_gzreadline(hdrstr, gzin)) > 0, "Empty file: %s", path);
  load_assert(hdrstr->buff[0] == '{', "Expected JSON header: %s", path);

  size_t i, prev_offset = 0;
  size_t num_curly_open = 0, num_curly_close = 0; // '{' '}'
  size_t num_brkt_open = 0, num_brkt_close = 0; // '[' ']'
  bool in_string = false, escape_char = false; // '\'

  while(1)
  {
    for(i = prev_offset; i < hdrstr->len; i++) {
      if(in_string) {
        if(escape_char) escape_char = false;
        else if(hdrstr->buff[i] == '\\') escape_char = true;
        else if(hdrstr->buff[i] == '"') in_string = false;
      }
      else if(hdrstr->buff[i] == '"') in_string = true;
      else if(hdrstr->buff[i] == '{') num_curly_open++;
      else if(hdrstr->buff[i] == '}') num_curly_close++;
      else if(hdrstr->buff[i] == '[') num_brkt_open++;
      else if(hdrstr->buff[i] == ']') num_brkt_close++;
    }
    prev_offset = hdrstr->len;

    if(num_curly_open == num_curly_close && num_brkt_open == num_brkt_close)
      break;

    // header is not finished yet
    load_assert(num_curly_open > num_curly_close, "'}' before '{': %s", path);
    load_assert(num_brkt_open >= num_brkt_close, "']' before '[': %s", path);
    load_assert(hdrstr->len < MAX_JSON_HDR_BYTES, "Large JSON header: %s", path);
    load_assert((nread = strbuf_gzreadline(hdrstr, gzin)) > 0,
                "Premature end of JSON header: %s", path);
  }
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

// Open file
// if cannot open file returns 0
// if fatal is true, exits on error
// if !fatal, returns -1 on error
// if successful creates a new GPathReader and returns 1
int gpath_reader_open2(GPathReader *file, char *path, const char *mode,
                       bool fatal)
{
  FileFilter *fltr = &file->fltr;

  if(!file_filter_open(fltr, path, mode, true, fatal)) return 0;

  // Load JSON header into file->hdrstr
  StrBuf *hdrstr = &file->hdrstr;
  if(hdrstr->buff == NULL) strbuf_alloc(hdrstr, 1024);
  _gpath_reader_load_json_hdr(fltr->gz, path, hdrstr);
  file->json = cJSON_Parse(hdrstr->buff);
  if(file->json == NULL) die("Invalid JSON header: %s", path);

  // DEV: validate json header
  _parse_json_header(file);

  // The following functions call die() if kmer_size or ncols header fields
  // are missing
  size_t kmer_size = gpath_reader_get_kmer_size(file);
  size_t filencols = _gpath_reader_get_filencols(file);
  file_filter_set_cols(fltr, filencols);

  // Check we can handle the kmer size
  db_graph_check_kmer_size(kmer_size, file->fltr.file_path.buff);

  return 1;
}


void gpath_reader_check(const GPathReader *file, size_t db_kmer_size,
                        size_t db_ncols)
{
  size_t kmer_size = gpath_reader_get_kmer_size(file);
  if(kmer_size != db_kmer_size) {
    die("Path file kmer size doesn't match [path: %s kmer: %zu, graph-kmer: %zu]",
        file->fltr.file_path.buff, kmer_size, db_kmer_size);
  }
  size_t used_ncols = file_filter_usedcols(&file->fltr);
  if(used_ncols > db_ncols) {
    die("Path file requires at least %zu colours [path: %s, curr colours: %zu]",
        used_ncols, file->fltr.orig_path.buff, db_ncols);
  }
}

int gpath_reader_open(GPathReader *file, char *path, bool fatal)
{
  return gpath_reader_open2(file, path, "r", fatal);
}

void gpath_reader_close(GPathReader *file)
{
  file_filter_close(&file->fltr);
  cJSON_Delete(file->json);
  strbuf_dealloc(&file->hdrstr);
  ctx_free(file->colours_json);
  memset(file, 0, sizeof(GPathReader));
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

  for(i = 0; i < line->len && char_is_acgt(line->buff[i]); i++) {}
  load_assert(i == kmer_size, "Bad kmer line: %s\n%s\n", path, line->buff);

  BinaryKmer bkmer, bkey;
  hkey_t hkey;
  bool found = false;

  bkmer = binary_kmer_from_str(line->buff, kmer_size);
  // check bkmer is kmer key
  bkey = binary_kmer_get_key(bkmer, kmer_size);
  load_assert(binary_kmers_are_equal(bkmer, bkey), "Bkmer not bkey: %s", path);

  if(dont_add_kmers) {
    hkey = hash_table_find(&db_graph->ht, bkey);
    load_assert(hkey != HASH_NOT_FOUND, "BKmer not already loaded: %s", path);
  } else {
    hkey = hash_table_find_or_insert(&db_graph->ht, bkey, &found);
  }

  // Parse number of paths
  numpstr = line->buff+kmer_size;
  load_assert(numpstr[0] == ' ', "Bad kmer line: %s", path);
  numpstr++;
  num_paths_exp = strtoul(numpstr, &endpstr, 10);
  load_assert(endpstr > numpstr && (*endpstr == '\0' || *endpstr == ' '),
              "Bad kmer line: %s\n%s\n", path, line->buff);

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

  orient = (line->buff[0] == 'F') ? FORWARD : REVERSE;
  load_assert(line->buff[1] == ' ', "Bad path line: %s", path);
  pstr = line->buff+2;

  // [nkmers]
  num_kmers = strtoul(pstr, &endpstr, 10);
  load_assert(endpstr > pstr && *endpstr == ' ', "Bad path line: %s", path);
  pstr = endpstr+1;

  // [njuncs]
  num_juncs = strtoul(pstr, &endpstr, 10);
  load_assert(endpstr > pstr && *endpstr == ' ', "Bad path line: %s", path);
  pstr = endpstr+1;

  // [nseen,nseen,...]
  endpstr = strchr(pstr, ' ');
  load_assert(endpstr != NULL && endpstr > pstr, "Bad path line: %s", path);
  
  for(i = 0; i < file->fltr.filencols; i++, pstr = endpstr+1)
  {
    size_t num_seen = strtoul(pstr, &endpstr, 10);
    nseenbuf[i] = MIN2(UINT8_MAX, num_seen);
    load_assert((i+1 < file->fltr.filencols && *endpstr == ',') || *endpstr == ' ',
                "Bad path line: %s\n%s", path, line->buff);
  }

  // [seq]
  endpstr = pstr;
  do {
    load_assert(char_is_acgt(*endpstr), "Bad path line: %s", path);
    endpstr++;
  } while(*endpstr && *endpstr != ' ');
  load_assert((size_t)(endpstr - pstr) == num_juncs, "Bad path line: %s", path);
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

  load_assert(gpset->entries.len == num_paths_exp,
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
  gzFile gzin = file->fltr.gz;
  const char *path = file->fltr.file_path.buff;

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

  int num_kmers_exp = _json_demand_int(file->json, "num_kmers",
                                       file->fltr.file_path.buff);
  if(num_kmers_exp < 0) die("Invlaid num_kmers: %i", num_kmers_exp);

  while(strbuf_reset_gzreadline(&line, gzin) > 0) {
    strbuf_chomp(&line);
    if(line.len > 0 && line.buff[0] != '#') {
      if(line.buff[0] == 'F' || line.buff[0] == 'R') {
        // Assume path line
        // status("Path line: %s", line.buff);
        load_assert(hkey != HASH_NOT_FOUND, "Path before kmer: %s", path);

        _gpath_reader_load_path_line(file, path, &line, nseenbuf, &seqbuf,
                                     &gpset);

        load_assert(gpset.entries.len <= num_paths_exp, "Too many paths: %s", path);
      }
      else {
        // Assume kmer line
        // status("Kmer line: %s", line.buff);

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

  load_assert(num_kmers == (size_t)num_kmers_exp,
              "num_kmers don't match (exp %zu vs %zu)",
              num_kmers, (size_t)num_kmers_exp);

  gpath_subset_dealloc(&subset0);
  gpath_subset_dealloc(&subset1);

  gpath_set_dealloc(&gpset);
  strbuf_dealloc(&line);
  byte_buf_dealloc(&seqbuf);
  ctx_free(nseenbuf);
}

// Get max mem required to load ctp files
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
    size_t pbytes = gpath_reader_get_path_bytes(&files[i]);
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

    if(strcmp(gname->buff, "undefined") == 0) {
      strbuf_set(gname, pname);
    } else if(strcmp(gname->buff, pname) != 0) {
      die("Graph/path sample names do not match [%zu->%zu] '%s' vs '%s'",
          fromcol, intocol, gname->buff, pname);
    }
  }
}
