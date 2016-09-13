#include "global.h"
#include "gpath_save.h"
#include "gpath_checks.h"
#include "gpath_set.h"
#include "gpath_subset.h"
#include "binary_seq.h"
#include "util.h"
#include "json_hdr.h"

const char ctp_explanation_comment[] =
"# This file was generated with McCortex\n"
"#   written by Isaac Turner <turner.isaac@gmail.com>\n"
"#   url: "MCCORTEX_URL"\n"
"# \n"
"# Comment lines begin with a # and are ignored, but must come after the header\n"
"# Format is:\n"
"#   [kmer] [num_paths] ...(ignored)\n"
"#   [FR] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] [seq=... juncpos=... ...]\n"
"#\n"
"# Columns are separated by a single space.\n"
"# Columns 1-4 are required ([FR],..,[juncs]) everything after than is optional\n"
"\n";

// {
//   <JSON header>
// }
// <KMER> <num> .. (ignored)
// [FR] (nkmers) [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. (ignored)
//
// Note: (nkmers) is only in versions <= 3

static void _gpath_save_contig_hist2json(cJSON *json_hists,
                                         const size_t *arr_counts,
                                         size_t arr_len)
{
  cJSON *lens = cJSON_CreateArray();
  cJSON *cnts = cJSON_CreateArray();
  size_t i;

  for(i = 1; i < arr_len; i++) {
    if(arr_counts[i]) {
      cJSON_AddItemToArray(lens, cJSON_CreateInt(i));
      cJSON_AddItemToArray(cnts, cJSON_CreateInt(arr_counts[i]));
    }
  }

  cJSON *hist = cJSON_CreateObject();
  cJSON_AddItemToObject(hist, "lengths", lens);
  cJSON_AddItemToObject(hist, "counts",  cnts);

  cJSON_AddItemToArray(json_hists, hist);
}

/**
 * Generate a JSON header object for a .ctp file
 * @param path        path to output file
 * @param cmdstr      name of the command being run, to be used to add @cmdhdr
 * @param cmdhdr      JSON header to add under current command->@cmdstr
 *                    If cmdstr and cmdhdr are both NULL they are ignored
 * @param contig_hist histgram of read contig lengths
 * @param hist_len    length of array contig_hist
 */
cJSON* gpath_save_mkhdr(const char *path,
                        const char *cmdstr, cJSON *cmdhdr,
                        cJSON **hdrs, size_t nhdrs,
                        const ZeroSizeBuffer *contig_hists, size_t ncols,
                        const dBGraph *db_graph)
{
  ctx_assert(!cmdstr == !cmdhdr);

  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;

  // using json_hdr_make_std() assumes the following
  ctx_assert(gpset->ncols == db_graph->num_of_cols);

  // Construct cJSON
  cJSON *json = cJSON_CreateObject();

  cJSON_AddStringToObject(json, "file_format", "ctp");
  cJSON_AddNumberToObject(json, "format_version", CTP_FORMAT_VERSION);

  // Add standard cortex header info, including the command being run
  json_hdr_make_std(json, path, hdrs, nhdrs, db_graph);

  // Get first command (this one), and command specific extra info
  if(cmdstr) {
    cJSON *cmd = json_hdr_get_curr_cmd(json, path);
    cJSON_AddItemToObject(cmd, cmdstr, cmdhdr);
  }

  // Paths info
  cJSON *paths = cJSON_CreateObject();
  cJSON_AddItemToObject(json, "paths", paths);

  // Add command specific header fields
  cJSON_AddNumberToObject(paths, "num_kmers_with_paths", gpstore->num_kmers_with_paths);
  cJSON_AddNumberToObject(paths, "num_paths", gpstore->num_paths);
  cJSON_AddNumberToObject(paths, "path_bytes", gpstore->path_bytes);

  // Add size distribution
  cJSON *json_hists = cJSON_CreateArray();
  cJSON_AddItemToObject(paths, "contig_hists", json_hists);

  size_t i;
  for(i = 0; i < ncols; i++)
    _gpath_save_contig_hist2json(json_hists, contig_hists[i].b, contig_hists[i].len);

  return json;
}


static inline void _gpath_save_flush(gzFile gzout, StrBuf *sbuf,
                                     pthread_mutex_t *outlock)
{
  pthread_mutex_lock(outlock);
  gzwrite(gzout, sbuf->b, sbuf->end);
  pthread_mutex_unlock(outlock);
  strbuf_reset(sbuf);
}

/**
 * Print paths to a string buffer. Paths are sorted before being written.
 *
 * @param hkey    All paths associated with hkey are written to the buffer
 * @param sbuf    paths are written this string buffer
 * @param subset  is a temp variable that is reused each time
 * @param nbuf    temporary buffer, if not NULL, used to add seq=... to output
 * @param jposbuf temporary buffer, if not NULL, used to add juncpos=... to output
 */
void gpath_save_sbuf(hkey_t hkey, StrBuf *sbuf, GPathSubset *subset,
                     dBNodeBuffer *nbuf, SizeBuffer *jposbuf,
                     const dBGraph *db_graph)
{
  ctx_assert(db_graph->num_of_cols == 1 || nbuf == NULL);
  ctx_assert(db_graph->num_of_cols == 1 || jposbuf == NULL);

  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;
  const size_t ncols = gpstore->gpset.ncols;
  GPath *first_gpath = gpath_store_fetch(gpstore, hkey);
  const GPath *gpath;
  size_t i, j, col;

  // Load and sort paths for given kmer
  gpath_subset_reset(subset);
  gpath_subset_load_llist(subset, first_gpath);
  gpath_subset_sort(subset);

  if(subset->list.len == 0) return;

  // Print "<kmer> <npaths>"
  BinaryKmer bkmer = db_graph->ht.table[hkey];
  char bkstr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(bkmer, db_graph->kmer_size, bkstr);

  // strbuf_sprintf(sbuf, "%s %zu\n", bkstr, subset->list.len);
  strbuf_append_strn(sbuf, bkstr, db_graph->kmer_size);
  strbuf_append_char(sbuf, ' ');
  strbuf_append_ulong(sbuf, subset->list.len);
  strbuf_append_char(sbuf, '\n');

  char orchar[2] = {0};
  orchar[FORWARD] = 'F';
  orchar[REVERSE] = 'R';
  const uint8_t *nseenptr;

  for(i = 0; i < subset->list.len; i++)
  {
    gpath = subset->list.b[i];
    nseenptr = gpath_set_get_nseen(gpset, gpath);

    // strbuf_sprintf(sbuf, "%c %zu %u %u", orchar[gpath->orient], klen,
    //                                      gpath->num_juncs, (uint32_t)nseenptr[0]);

    strbuf_append_char(sbuf, orchar[gpath->orient]);
    strbuf_append_char(sbuf, ' ');
    strbuf_append_ulong(sbuf, gpath->num_juncs);
    strbuf_append_char(sbuf, ' ');
    strbuf_append_ulong(sbuf, nseenptr[0]);

    for(col = 1; col < ncols; col++) {
      // strbuf_sprintf(sbuf, ",%u", (uint32_t)nseenptr[col]);
      strbuf_append_char(sbuf, ',');
      strbuf_append_ulong(sbuf, nseenptr[col]);
    }

    strbuf_append_char(sbuf, ' ');
    strbuf_ensure_capacity(sbuf, sbuf->end + gpath->num_juncs + 2);
    binary_seq_to_str(gpath->seq, gpath->num_juncs, sbuf->b+sbuf->end);
    sbuf->end += gpath->num_juncs;

    if(nbuf)
    {
      // Trace this path through the graph
      // First, find a colour this path is in
      for(col = 0; col < ncols && !gpath_has_colour(gpath, ncols, col); col++) {}
      if(col == ncols) die("path is not in any colours");

      dBNode node = {.key = hkey, .orient = gpath->orient};
      db_node_buf_reset(nbuf);
      if(jposbuf) size_buf_reset(jposbuf); // indices of junctions in nbuf
      gpath_fetch(node, gpath, nbuf, jposbuf, col, db_graph);

      strbuf_append_str(sbuf, " seq=");
      strbuf_ensure_capacity(sbuf, sbuf->end + db_graph->kmer_size + nbuf->len);
      sbuf->end += db_nodes_to_str(nbuf->b, nbuf->len, db_graph,
                                   sbuf->b+sbuf->end);

      if(jposbuf) {
        strbuf_append_str(sbuf, " juncpos=");
        strbuf_append_ulong(sbuf, jposbuf->b[0]);

        for(j = 1; j < jposbuf->len; j++) {
          strbuf_append_char(sbuf, ',');
          strbuf_append_ulong(sbuf, jposbuf->b[j]);
        }
      }
    }

    strbuf_append_char(sbuf, '\n');
  }
}

// @subset is a temp variable that is reused each time
// @sbuf   is a temp variable that is reused each time
static inline int _gpath_gzsave_node(hkey_t hkey,
                                     StrBuf *sbuf, GPathSubset *subset,
                                     dBNodeBuffer *nbuf, SizeBuffer *jposbuf,
                                     gzFile gzout, pthread_mutex_t *outlock,
                                     const dBGraph *db_graph)
{
  gpath_save_sbuf(hkey, sbuf, subset, nbuf, jposbuf, db_graph);

  if(sbuf->end > DEFAULT_IO_BUFSIZE)
    _gpath_save_flush(gzout, sbuf, outlock);

  return 0; // => keep iterating
}

typedef struct
{
  size_t nthreads;
  bool save_seq; // write seq=... juncpos=...
  gzFile gzout;
  pthread_mutex_t *outlock;
  dBGraph *db_graph;
} GPathSaving;

static void gpath_save_thread(void *arg, size_t threadid)
{
  GPathSaving *save = (GPathSaving*)arg;
  const dBGraph *db_graph = save->db_graph;

  GPathSubset subset;
  StrBuf sbuf;

  gpath_subset_alloc(&subset);
  gpath_subset_init(&subset, &save->db_graph->gpstore.gpset);
  strbuf_alloc(&sbuf, 2 * DEFAULT_IO_BUFSIZE);

  dBNodeBuffer nbuf;
  SizeBuffer jposbuf;
  db_node_buf_alloc(&nbuf, 1024);
  size_buf_alloc(&jposbuf, 256);

  HASH_ITERATE_PART(&db_graph->ht, threadid, save->nthreads,
                    _gpath_gzsave_node,
                    &sbuf, &subset,
                    save->save_seq ? &nbuf : NULL, save->save_seq ? &jposbuf : NULL,
                    save->gzout, save->outlock,
                    db_graph);

  _gpath_save_flush(save->gzout, &sbuf, save->outlock);

  db_node_buf_dealloc(&nbuf);
  size_buf_dealloc(&jposbuf);
  gpath_subset_dealloc(&subset);
  strbuf_dealloc(&sbuf);
}

/**
 * Save paths to a file.
 * @param gzout         gzFile to write to
 * @param path          path of output file
 * @param save_path_seq if true, save seq= and juncpos= for links, requires
 *                      exactly one colour in the graph
 * @param hdrs is array of JSON headers of input files
 */
void gpath_save(gzFile gzout, const char *path,
                size_t nthreads, bool save_path_seq,
                const char *cmdstr, cJSON *cmdhdr,
                cJSON **hdrs, size_t nhdrs,
                const ZeroSizeBuffer *contig_hists, size_t ncols,
                dBGraph *db_graph)
{
  ctx_assert(nthreads > 0);
  ctx_assert(gpath_set_has_nseen(&db_graph->gpstore.gpset));
  ctx_assert(ncols == db_graph->gpstore.gpset.ncols);
  ctx_assert(!save_path_seq || db_graph->num_of_cols == 1); // save_path => 1 colour

  char npaths_str[50];
  ulong_to_str(db_graph->gpstore.num_paths, npaths_str);

  status("Saving %s paths to: %s", npaths_str, path);
  status("  using %zu threads", nthreads);

  // Write header
  cJSON *json = gpath_save_mkhdr(path, cmdstr, cmdhdr, hdrs, nhdrs,
                                 contig_hists, ncols, db_graph);
  json_hdr_gzprint(json, gzout);
  cJSON_Delete(json);

  // Print comments about the format
  gzputs(gzout, ctp_explanation_comment);

  // Multithreaded
  pthread_mutex_t outlock;
  if(pthread_mutex_init(&outlock, NULL) != 0) die("Mutex init failed");

  GPathSaving save = {.nthreads = nthreads,
                      .save_seq = save_path_seq,
                      .gzout = gzout,
                      .outlock = &outlock,
                      .db_graph = db_graph};

  // Iterate over kmers writing paths
  util_multi_thread(&save, nthreads, gpath_save_thread);
  pthread_mutex_destroy(&outlock);
  status("[GPathSave] Graph paths saved to %s", path);
}
