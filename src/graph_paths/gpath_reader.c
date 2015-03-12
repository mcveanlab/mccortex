#include "global.h"
#include "gpath_reader.h"
#include "file_util.h"
#include "util.h"
#include "hash_mem.h"
#include "common_buffers.h"
#include "str_parsing.h" // comma_list_to_array()
#include "binary_seq.h"
#include "binary_kmer.h"
#include "db_node.h"
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

size_t gpath_reader_get_kmer_size(const GPathReader *file)
{
  return json_hdr_get_kmer_size(file->json, file->fltr.path.b);
}

size_t gpath_reader_get_num_kmers(const GPathReader *file)
{
  cJSON *paths = json_hdr_get_paths(file->json, file->fltr.path.b);
  return json_hdr_demand_uint(paths, "num_kmers_with_paths", file->fltr.path.b);
}

size_t gpath_reader_get_num_paths(const GPathReader *file)
{
  cJSON *paths = json_hdr_get_paths(file->json, file->fltr.path.b);
  return json_hdr_demand_uint(paths, "num_paths", file->fltr.path.b);
}

size_t gpath_reader_get_path_bytes(const GPathReader *file)
{
  cJSON *paths = json_hdr_get_paths(file->json, file->fltr.path.b);
  return json_hdr_demand_uint(paths, "path_bytes", file->fltr.path.b);
}

static size_t _gpath_reader_get_filencols(const GPathReader *file)
{
  return json_hdr_get_ncols(file->json, file->fltr.path.b);
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

void gpath_reader_load_contig_hist(cJSON *json_root, const char *path,
                                   size_t fromcol, ZeroSizeBuffer *hist)
{
  size_t i;
  cJSON *json_fmt = json_hdr_get(json_root, "file_format", cJSON_String, path);
  cJSON *hdr_paths, *contig_hists, *json_hist, *json_lens, *json_cnts;

  if(strcmp(json_fmt->valuestring,"ctp") != 0) {
    die("File is not CTP '%s' [%s]", json_fmt->valuestring, path);
  }

  // File is ctp
  hdr_paths = json_hdr_get(json_root, "paths", cJSON_Object, path);
  contig_hists = json_hdr_get(hdr_paths, "contig_hists", cJSON_Array, path);

  json_hist = cJSON_GetArrayItem(contig_hists, fromcol);
  if(json_hist == NULL) die("No hist for colour %zu [path: %s]", fromcol, path);

  json_lens = json_hdr_get(json_hist, "lengths", cJSON_Array, path);
  json_cnts = json_hdr_get(json_hist, "counts" , cJSON_Array, path);

  // Loop over entries
  cJSON *len = json_lens->child, *cnt = json_cnts->child;
  for(i = 0; len && cnt; i++, len = len->next, cnt = cnt->next) {
    if(len->type != cJSON_Number || cnt->type != cJSON_Number ||
       len->valueint < 0 || cnt->valueint < 0) {
      die("JSON array entry is not a number");
    }
    ctx_assert2(len->valueint < 1<<30, "%li is too big", len->valueint);
    zsize_buf_resize(hist, len->valueint+1);
    hist->b[len->valueint] += cnt->valueint;
  }

  if(len || cnt) {
    die("Contig histogram arrays not the same lengths");
  }
}

void gpath_reader_get_max_contig_lens(GPathReader *file, size_t *max_contigs)
{
  cJSON *hdr_paths, *contig_hists, *json_hist, *json_lens, *json_cnts;
  size_t i, max_contig, fromcol, intocol, ncols = file_filter_num(&file->fltr);
  const char *path = file_filter_path(&file->fltr);

  hdr_paths = json_hdr_get(file->json, "paths", cJSON_Object, path);
  contig_hists = json_hdr_get(hdr_paths, "contig_hists", cJSON_Array, path);

  for(i = 0; i < ncols; i++)
  {
    fromcol = file_filter_fromcol(&file->fltr, i);
    intocol = file_filter_intocol(&file->fltr, i);
    json_hist = cJSON_GetArrayItem(contig_hists, fromcol);
    if(json_hist == NULL) die("No hist for colour %zu [path: %s]", fromcol, path);

    json_lens = json_hdr_get(json_hist, "lengths", cJSON_Array, path);
    json_cnts = json_hdr_get(json_hist, "counts" , cJSON_Array, path);

    // Loop over entries
    cJSON *len = json_lens->child, *cnt = json_cnts->child;
    for(max_contig = 0; len && cnt; len = len->next, cnt = cnt->next) {
      if(len->type != cJSON_Number || cnt->type != cJSON_Number ||
         len->valueint < 0 || cnt->valueint < 0) {
        die("JSON array entry is not a number");
      }
      if(cnt->valueint > 0) max_contig = MAX2(max_contig, (size_t)len->valueint);
    }

    if(len || cnt) die("Contig histogram arrays not the same lengths");

    max_contigs[intocol] = MAX2(max_contigs[intocol], max_contig);
  }
}

/**
 * Parse values from file->json into fields in file
 */
void _parse_json_header(GPathReader *file)
{
  cJSON *graph, *colours, *sample;

  graph = cJSON_GetObjectItem(file->json, "graph");
  if(graph == NULL || graph->type != cJSON_Object)
    die("No 'graph' array entry in header");

  colours = cJSON_GetObjectItem(graph, "colours");
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
void gpath_reader_open2(GPathReader *file, const char *path, const char *mode,
                        size_t into_offset)
{
  FileFilter *fltr = &file->fltr;
  file_filter_open(fltr, path); // calls die() on error

  file->gz = futil_gzopen(fltr->path.b, mode);
  strm_buf_alloc(&file->strmbuf, 4*ONE_MEGABYTE);

  // Temporary variable for loading
  strbuf_alloc(&file->line, 1024);
  size_buf_alloc(&file->numbuf, 16);

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
  file_filter_set_cols(fltr, filencols, into_offset);

  // Check we can handle the kmer size
  db_graph_check_kmer_size(kmer_size, file->fltr.path.b);
}

void gpath_reader_open(GPathReader *file, const char *path)
{
  gpath_reader_open2(file, path, "r", 0);
}

void gpath_reader_close(GPathReader *file)
{
  if(file->gz) gzclose(file->gz);
  strm_buf_dealloc(&file->strmbuf);
  strbuf_dealloc(&file->line);
  size_buf_dealloc(&file->numbuf);
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

// Reads line <kmer> <num_links>
// Calls die() on error
// Returns true unless end of file
bool gpath_reader_read_kmer(GPathReader *file, StrBuf *kmer, size_t *num_links)
{
  strbuf_reset(kmer);
  *num_links = 0;

  const char *path = file_filter_path(&file->fltr);
  int c;
  char *space;

  while((c = gzgetc_buf(file->gz, &file->strmbuf)) != -1)
  {
    if(c == '#') gzskipline_buf(file->gz, &file->strmbuf);
    else if(c != '\n') {
      strbuf_append_char(kmer, c);
      strbuf_gzreadline_buf(kmer, file->gz, &file->strmbuf);
      strbuf_chomp(kmer);
      if(!char_is_acgt(c) ||
         (space = strchr(kmer->b, ' ')) == NULL ||
         !parse_entire_size(space+1, num_links))
      {
        die("Bad kmer line [%s]: %s", path, kmer->b);
      }
      strbuf_resize(kmer, space - kmer->b);
      return true;
    }
  }

  return false;
}

#define bad_link_line(path,line) die("Bad link line [%s]: %s", path, (line)->b);

// [FR] [nkmers] [njuncs] [nseen0,nseen1,...] [juncs:ACAGT] ([seq=] [juncpos=])?
void parse_link_line(GPathReader *file, const StrBuf *line, SizeBuffer *numbuf,
                     bool *fw, size_t *kdist, size_t *njuncs,
                     SizeBuffer *counts, StrBuf *juncs,
                     StrBuf *seq, SizeBuffer *juncpos)
{
  const char *path = file_filter_path(&file->fltr);
  size_t i, fromcol, intocol;
  char *end = NULL;

  // First first 5 required columns
  const char *cols[5] = {NULL};
  for(cols[0] = line->b, i = 1; i < 5; i++) {
    if((cols[i] = strchr(cols[i-1], ' ')) == NULL) bad_link_line(path,line);
    cols[i]++;
  }
  // Get column lengths
  size_t collens[5] = {0};
  for(i = 0; i < 4; i++) collens[i] = cols[i+1]-cols[i]-1;
  collens[4] = strendc(cols[4], ' ') - cols[4];

  // 0:[FR]
  char c = line->b[0];
  if((c != 'F' && c != 'R') || collens[0] != 1) bad_link_line(path,line);
  *fw = (c == 'F');

  // 1:[nkmers]
  *kdist = strtoul(cols[1], &end, 10);
  if(end != cols[1]+collens[1]) bad_link_line(path,line);

  // 2:[njuncs]
  *njuncs = strtoul(cols[2], &end, 10);
  if(end != cols[2]+collens[2]) bad_link_line(path,line);

  // 3:[nseen0,nseen1,...]
  size_buf_reset(numbuf);
  size_buf_reset(counts);
  if(comma_list_to_array(cols[3], numbuf) != (int)collens[3])
    bad_link_line(path,line);

  // Use filter - append zeros first
  size_buf_push_zero(counts, file_filter_into_ncols(&file->fltr));
  for(i = 0; i < numbuf->len; i++) {
    fromcol = file_filter_fromcol(&file->fltr, i);
    intocol = file_filter_intocol(&file->fltr, i);
    counts->b[intocol] += numbuf->b[fromcol];
  }

  // 4:[juncs:ACAGA]
  strbuf_reset(juncs);
  strbuf_append_strn(juncs, cols[4], collens[4]);

  // Parse optional tags
  const char *txt;
  size_t txtlen;

  // [seq=ACACA]
  strbuf_reset(seq);
  if(str_find_tag(line->b, "seq=", &txt, &txtlen))
  {
    strbuf_append_strn(seq, txt, txtlen);
    for(i = 0; i < txtlen; i++)
      if(!char_is_acgt(txt[i]))
        die("Cannot parse seq= [%s]: %s", path, line->b);
  }

  // [juncpos=12,23]
  size_buf_reset(juncpos);
  if(str_find_tag(line->b, "juncpos=", &txt, &txtlen) &&
     comma_list_to_array(txt, juncpos) != (int)txtlen)
  {
    die("Cannot parse juncpos= [%s]: %s", path, line->b);
  }

  //
  // Some sanity checks
  //
  if(juncs->end != *njuncs) die("Differing lengths: %s", line->b);
  if(*kdist < *njuncs+2) die("kdist too short");

  if(juncpos->len > 0) {
    size_t last = juncpos->b[juncpos->len-1];
    if(juncpos->len != *njuncs) die("Mismatch %zu vs %zu", juncpos->len, *njuncs);
    if(last+2 != *kdist) die("Mismatch %zu vs %zu", last+2, *kdist);
  }

  if(juncpos->len && seq->end) {
    if(juncpos->b[juncpos->len-1] >= seq->end) die("Seq too short");
    size_t p, ksize = seq->end - (juncpos->b[juncpos->len-1] + 1);
    if(!(ksize & 1)) die("kmer_size not odd");
    for(i = 0; i < juncpos->len; i++) {
      p = juncpos->b[i];
      if(p+ksize > seq->end || seq->b[ksize+p] != juncs->b[i])
        die("Bad entry");
    }
  }
}

// Reads line [FR] <num_links>
// Calls die() on error
// Returns true unless end of link entries
bool gpath_reader_read_link(GPathReader *file,
                            bool *fw, size_t *kdist, size_t *njuncs,
                            SizeBuffer *countbuf, StrBuf *juncs,
                            StrBuf *seq, SizeBuffer *juncpos)
{
  int c;
  StrBuf *line = &file->line;
  strbuf_reset(line);

  while((c = gzgetc_buf(file->gz, &file->strmbuf)) != -1)
  {
    if(char_is_acgt(c)) return false; // Hit kmer line
    else if(c == '#') gzskipline_buf(file->gz, &file->strmbuf);
    else if(c != '\n') {
      strbuf_append_char(line, c);
      strbuf_gzreadline_buf(line, file->gz, &file->strmbuf);
      strbuf_chomp(line);
      parse_link_line(file, line, &file->numbuf,
                      fw, kdist, njuncs, countbuf, juncs,
                      seq, juncpos);
      return true;
    }
  }

  return false;
}

typedef struct
{
  BinaryKmer bkey;
  hkey_t hkey;
  int kmer_flags; // how to load this kmer
  size_t num_paths_exp, num_paths_seen;
  bool found; // bkey was already in the hash table
  bool looked_up; // Whether we have tried to find the key in the graph
} LoadPathKmer;

static void _validate_new_path(LoadPathKmer *load, const char *path,
                               uint8_t *nseen, size_t nseencols,
                               dBGraph *db_graph)
{
  if(load->looked_up) return;

  bool found = false;

  switch(load->kmer_flags) {
    case GPATH_ADD_MISSING_KMERS:
      load->hkey = hash_table_find_or_insert(&db_graph->ht, load->bkey, &found);
      break;
    case GPATH_DIE_MISSING_KMERS:
      load->hkey = hash_table_find(&db_graph->ht, load->bkey);
      if(load->hkey == HASH_NOT_FOUND) {
        char bkmerstr[MAX_KMER_SIZE+1];
        binary_kmer_to_str(load->bkey, db_graph->kmer_size, bkmerstr);
        die("BKmer not already loaded: %s [%s]", bkmerstr, path);
      }
      found = true;
      break;
    case GPATH_SKIP_MISSING_KMERS:
      load->hkey = hash_table_find(&db_graph->ht, load->bkey);
      found = (load->hkey != HASH_NOT_FOUND);
      break;
    default: die("Bad switch value: %i", load->kmer_flags);
  }

  // If inserted, add to colour
  size_t i;
  if(!found && load->hkey != HASH_NOT_FOUND && db_graph->node_in_cols) {
    for(i = 0; i < nseencols; i++)
      db_node_or_col(db_graph, load->hkey, i, nseen[i] > 0);
  }

  load->found = found;
  load->looked_up = true;
}

/**
 * <KMER> <num> .. (ignored)
 * @param kmer_flags must be one of:
 *   * GPATH_ADD_MISSING_KMERS - add kmers to the graph before loading path
 *   * GPATH_DIE_MISSING_KMERS - die with error if cannot find kmer
 *   * GPATH_SKIP_MISSING_KMERS - skip paths where kmer is not in graph
 */
static void _gpath_reader_load_kmer_line(const char *path,
                                         StrBuf *line, int kmer_flags,
                                         LoadPathKmer *result,
                                         dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  char *numpstr, *endpstr;
  size_t i, num_paths_exp = 0;
  BinaryKmer bkmer, bkey;

  for(i = 0; i < line->end && char_is_acgt(line->b[i]); i++) {}
  load_check(i == kmer_size, "Bad kmer line: %s\n'%s'\n", path, line->b);

  bkmer = binary_kmer_from_str(line->b, kmer_size);
  // check bkmer is kmer key
  bkey = binary_kmer_get_key(bkmer, kmer_size);
  load_check(binary_kmers_are_equal(bkmer, bkey), "Bkmer not bkey: %s", path);

  // Parse number of paths
  numpstr = line->b+kmer_size;
  load_check(numpstr[0] == ' ', "Bad kmer line: %s", path);
  numpstr++;
  num_paths_exp = strtoul(numpstr, &endpstr, 10);
  load_check(endpstr > numpstr && (*endpstr == '\0' || *endpstr == ' '),
              "Bad kmer line: %s\n%s\n", path, line->b);

  LoadPathKmer load = {.bkey = bkey, .hkey = HASH_NOT_FOUND,
                       .found = false, .looked_up = false,
                       .kmer_flags = kmer_flags,
                       .num_paths_exp = num_paths_exp,
                       .num_paths_seen = 0};

  memcpy(result, &load, sizeof(LoadPathKmer));
}

/**
 * Format: [FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. (ignored)
 * @param line        buffer holding input line to parse
 * @param tmp_nseen1  temporary memory of length file->fltr.filencols
 * @param tmp_nseen2  temporary memory of length intoncols(file->fltr)
 * @param tmp_seqbuf  temporary memory buffer
 * @param gpset       GPathSet to add new path to
 * @param db_graph    We look up and add kmers to this graph
 */
static void _gpath_reader_load_path_line(GPathReader *file, StrBuf *line,
                                         uint8_t *tmp_nseen1,
                                         uint8_t *tmp_nseen2,
                                         ByteBuffer *tmp_seqbuf,
                                         LoadPathKmer *load_kmer,
                                         size_t *max_contigs,
                                         GPathSet *gpset,
                                         dBGraph *db_graph)
{
  const char *path = file_filter_path(&file->fltr);
  size_t into_ncols = file_filter_into_ncols(&file->fltr);
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

  load_check(num_kmers > num_juncs, "%zu %zu", num_kmers, num_juncs);

  // [nseen,nseen,...]
  endpstr = strchr(pstr, ' ');
  load_check(endpstr != NULL && endpstr > pstr, "Bad path line: %s", path);

  for(i = 0; i < file->fltr.filencols; i++, pstr = endpstr+1)
  {
    size_t num_seen = strtoul(pstr, &endpstr, 10);
    tmp_nseen1[i] = MIN2(UINT8_MAX, num_seen);
    load_check((i+1 < file->fltr.filencols && *endpstr == ',') || *endpstr == ' ',
               "Bad path line: %s\n%s", path, line->b);
  }

  load_check(num_kmers < GPATH_MAX_KMERS, "%zu", (size_t)num_kmers);
  load_check(num_juncs < GPATH_MAX_JUNCS, "%zu", (size_t)num_juncs);

  // [seq]
  endpstr = pstr;
  do {
    load_check(char_is_acgt(*endpstr), "Bad path line: %s", path);
    endpstr++;
  } while(*endpstr && *endpstr != ' ');
  load_check((size_t)(endpstr - pstr) == num_juncs, "Bad path line: %s", path);
  byte_buf_capacity(tmp_seqbuf, (num_juncs+3)/4);
  binary_seq_from_str(pstr, num_juncs, tmp_seqbuf->b);

  // Filter colours
  bool path_in_cols = false;
  size_t contig_len = num_kmers + db_graph->kmer_size - 1;
  size_t fromcol, intocol;
  memset(tmp_nseen2, 0, sizeof(uint8_t) * into_ncols);

  for(i = 0; i < file_filter_num(&file->fltr); i++)
  {
    fromcol = file_filter_fromcol(&file->fltr, i);
    intocol = file_filter_intocol(&file->fltr, i);
    tmp_nseen2[intocol] = MIN2((size_t)UINT8_MAX,
                               (size_t)tmp_nseen2[intocol] + tmp_nseen1[fromcol]);

    if(tmp_nseen2[intocol] > 0) {
      path_in_cols = true;
      if(contig_len > max_contigs[intocol]) die("Invalid CTP contig histogram");
    }
  }

  if(path_in_cols)
  {
    // Load kmer
    _validate_new_path(load_kmer, path, tmp_nseen2, into_ncols, db_graph);

    if(load_kmer->hkey != HASH_NOT_FOUND)
    {
      // Add to GPathSet
      GPathNew newgpath = {.seq = tmp_seqbuf->b,
                           .colset = NULL, .nseen = NULL,
                           .orient = orient,
                           .klen = num_kmers,
                           .num_juncs = num_juncs};

      GPath *gpath = gpath_set_add_mt(gpset, newgpath);

      ctx_assert(into_ncols <= gpset->ncols);

      // Update nseen / colour bitset
      uint8_t *nseen = gpath_set_get_nseen(gpset, gpath);
      uint8_t *colset = gpath_get_colset(gpath, gpset->ncols);

      for(i = 0; i < into_ncols; i++) {
        nseen[i] = MIN2((size_t)UINT8_MAX, (size_t)nseen[i] + tmp_nseen2[i]);
        bitset_or(colset, i, nseen[i] > 0);
      }
    }
  }

  load_kmer->num_paths_seen++;
}

// @subset0 and @subset1 are temporary memory to be used in the loading
// of paths. Both are reset before being used.
// @return number of paths added
static size_t _load_paths_from_set(dBGraph *db_graph, GPathSet *gpset,
                                   GPathSubset *subset0, GPathSubset *subset1,
                                   const LoadPathKmer *load_kmer,
                                   const char *file_path)
{
  GPathStore *gpstore = &db_graph->gpstore;
  GPathHash *gphash = &db_graph->gphash;

  load_check(load_kmer->num_paths_exp == load_kmer->num_paths_seen,
             "Too many/few paths: %s (exp %zu vs act %zu)",
             file_path, load_kmer->num_paths_exp, load_kmer->num_paths_seen);

  size_t i;
  hkey_t hkey = load_kmer->hkey;

  gpath_subset_init(subset0, &gpstore->gpset);
  gpath_subset_init(subset1, gpset);

  GPath *kmer_paths = gpath_store_fetch(gpstore, hkey);
  gpath_subset_load_llist(subset0, kmer_paths);
  gpath_subset_load_set(subset1);
  gpath_subset_rmdup(subset1);

  // Merge entries from subset1 into subset0
  // false => do not remove duplicates
  gpath_subset_merge(subset0, subset1, false);

  bool found;

  // Copy remaining entries across
  GPathNew newgp;
  if(db_graph_has_path_hash(db_graph)) {
    for(i = 0; i < subset1->list.len; i++) {
      newgp = gpath_set_get(gpset, subset1->list.b[i]);
      gpath_hash_find_or_insert_mt(gphash, hkey, newgp, &found);
    }
  } else {
    for(i = 0; i < subset1->list.len; i++) {
      newgp = gpath_set_get(gpset, subset1->list.b[i]);
      gpath_store_add_mt(gpstore, hkey, newgp);
    }
  }

  gpath_set_reset(gpset);

  return subset1->list.len;
}

/**
 * @param kmer_flags must be one of:
 *   * GPATH_ADD_MISSING_KMERS - add kmers to the graph before loading path
 *   * GPATH_DIE_MISSING_KMERS - die with error if cannot find kmer
 *   * GPATH_SKIP_MISSING_KMERS - skip paths where kmer is not in graph
 */
void gpath_reader_load(GPathReader *file, int kmer_flags, dBGraph *db_graph)
{
  const char *path = file_filter_path(&file->fltr);

  file_filter_status(&file->fltr);

  // Read kmer lines
  // <KMER> <num>
  // [FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. ignored
  StrBuf line;
  strbuf_alloc(&line, 2048);

  size_t into_ncols = file_filter_into_ncols(&file->fltr);
  uint8_t *nseenbuf1 = ctx_calloc(file->fltr.filencols, sizeof(uint8_t));
  uint8_t *nseenbuf2 = ctx_calloc(into_ncols, sizeof(uint8_t));

  size_t *max_contigs = ctx_calloc(into_ncols, sizeof(size_t));
  gpath_reader_get_max_contig_lens(file, max_contigs);

  ByteBuffer seqbuf;
  byte_buf_alloc(&seqbuf, 64);

  // Load paths into this temporary set for each kmer
  GPathSet gpset;
  gpath_set_alloc(&gpset, db_graph->num_of_cols, ONE_MEGABYTE, true, true);

  GPathSubset subset0, subset1;
  gpath_subset_alloc(&subset0);
  gpath_subset_alloc(&subset1);

  LoadPathKmer load_kmer;
  size_t num_kmers = 0, num_kmers_exp = gpath_reader_get_num_kmers(file);
  size_t num_paths_loaded = 0, num_kmers_loaded = 0;

  while(strbuf_gzreadline_buf(&line, file->gz, &file->strmbuf) > 0) {
    strbuf_chomp(&line);
    if(line.end > 0 && line.b[0] != '#') {
      if(line.b[0] == 'F' || line.b[0] == 'R') {
        // Assume path line
        // status("Path line: %s", line.b);
        load_check(num_kmers != 0, "Path before kmer: %s", path);

        _gpath_reader_load_path_line(file, &line,
                                     nseenbuf1, nseenbuf2, &seqbuf,
                                     &load_kmer, max_contigs,
                                     &gpset, db_graph);
      }
      else {
        // Assume kmer line
        // Load previous paths if this is not the first kmer
        if(num_kmers > 0 && load_kmer.hkey != HASH_NOT_FOUND) {
          num_kmers_loaded += (gpset.entries.len > 0);
          num_paths_loaded += _load_paths_from_set(db_graph, &gpset,
                                                   &subset0, &subset1,
                                                   &load_kmer, path);
        }

        _gpath_reader_load_kmer_line(path, &line, kmer_flags, &load_kmer, db_graph);

        num_kmers++;
      }
    }
    strbuf_reset(&line);
  }

  if(num_kmers > 0 && load_kmer.hkey != HASH_NOT_FOUND) {
    num_kmers_loaded += (gpset.entries.len > 0);
    num_paths_loaded += _load_paths_from_set(db_graph, &gpset,
                                             &subset0, &subset1,
                                             &load_kmer, path);
  }

  load_check(num_kmers == (size_t)num_kmers_exp,
             "num_kmers don't match (exp %zu vs %zu)",
             num_kmers, (size_t)num_kmers_exp);

  // Print status update
  char npaths_str[50], nkmers_str[50];
  ulong_to_str(num_paths_loaded, npaths_str);
  ulong_to_str(num_kmers_loaded, nkmers_str);
  status("Loaded %s paths from %s kmers", npaths_str, nkmers_str);

  gpath_subset_dealloc(&subset0);
  gpath_subset_dealloc(&subset1);

  gpath_set_dealloc(&gpset);
  strbuf_dealloc(&line);
  byte_buf_dealloc(&seqbuf);
  ctx_free(nseenbuf1);
  ctx_free(nseenbuf2);
  ctx_free(max_contigs);
}

void gpath_reader_load_sample_names(const GPathReader *file, dBGraph *db_graph)
{
  const FileFilter *fltr = &file->fltr;

  size_t i, intocol, fromcol;
  StrBuf *gname;
  const char *pname;

  for(i = 0; i < file_filter_num(fltr); i++)
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
  size_t path_bytes, sum_path_bytes = 0;

  for(i = 0; i < nfiles; i++) {
    npaths = gpath_reader_get_num_paths(&files[i]);
    max_npaths = MAX2(max_npaths, npaths);
    sum_npaths += npaths;
    path_bytes = gpath_reader_get_path_bytes(&files[i]);
    sum_path_bytes += path_bytes;
  }

  npaths = sum_mem <= mem ? sum_npaths : 0;
  status("[GPathReader] need %zu paths %zu bytes", npaths, sum_path_bytes);

  gpath_store_alloc(&db_graph->gpstore,
                    db_graph->num_of_cols,
                    db_graph->ht.capacity,
                    npaths, mem,
                    count_nseen, false);
}
