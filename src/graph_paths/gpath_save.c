#include "global.h"
#include "gpath_save.h"
#include "gpath_set.h"
#include "gpath_subset.h"
#include "binary_seq.h"
#include "util.h"
#include "json_hdr.h"

// {
//   <JSON header>
// }
// <KMER> <num> .. (ignored)
// [FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. (ignored)

static void _gpath_save_contig_hist2json(cJSON *json_paths,
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

  cJSON_AddItemToObject(json_paths, "contig_lengths", lens);
  cJSON_AddItemToObject(json_paths, "contig_counts", cnts);
}

/**
 * @param contig_hist histgram of read contig lengths
 * @param hist_len    length of array contig_hist
 */
static void _gpath_save_hdr(gzFile gzout, const char *path,
                            cJSON **hdrs, size_t nhdrs,
                            const size_t *contig_hist, size_t hist_len,
                            const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;

  // Construct cJSON
  cJSON *json = cJSON_CreateObject();

  cJSON_AddStringToObject(json, "file_format", "ctp");
  cJSON_AddNumberToObject(json, "format_version", 3);

  // using json_hdr_add_std assumes the following
  ctx_assert(gpset->ncols == db_graph->num_of_cols);

  // Add standard cortex header info
  json_hdr_add_std(json, path, hdrs, nhdrs, db_graph);

  // Get first command (this one)
  // cJSON *cmd = json_hdr_get_curr_cmd(json);

  // Paths info
  cJSON *paths = cJSON_CreateObject();
  cJSON_AddItemToObject(json, "paths", paths);

  // Add command specific header fields
  cJSON_AddNumberToObject(paths, "num_kmers_with_paths", gpstore->num_kmers_with_paths);
  cJSON_AddNumberToObject(paths, "num_paths", gpstore->num_paths);
  cJSON_AddNumberToObject(paths, "path_bytes", gpstore->path_bytes);

  // Add size distribution
  // _gpath_save_contig_hist_merge(paths, hdrs, nhdrs, contig_hist, hist_len);
  _gpath_save_contig_hist2json(paths, contig_hist, hist_len);

  // Write header to file
  json_hdr_gzprint(json, gzout);

  cJSON_Delete(json);
}

static inline void _gpath_save_flush(gzFile gzout, StrBuf *sbuf,
                                     pthread_mutex_t *outlock)
{
  pthread_mutex_lock(outlock);
  gzwrite(gzout, sbuf->b, sbuf->end);
  pthread_mutex_unlock(outlock);
  strbuf_reset(sbuf);
}

// @subset is a temp variable that is reused each time
// @sbuf   is a temp variable that is reused each time
static inline int _gpath_save_node(hkey_t hkey,
                                   gzFile gzout, pthread_mutex_t *outlock,
                                   StrBuf *sbuf, GPathSubset *subset,
                                   const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;
  const size_t ncols = gpstore->gpset.ncols;
  GPath *first_gpath = gpath_store_fetch(gpstore, hkey);
  const GPath *gpath;
  size_t i, col;

  // Load and sort paths for given kmer
  gpath_subset_reset(subset);
  gpath_subset_load_llist(subset, first_gpath);
  gpath_subset_sort(subset);

  if(subset->list.len == 0) return 0; // => keep iterating

  // Print "<kmer> <npaths>"
  BinaryKmer bkmer = db_graph->ht.table[hkey];
  char bkstr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(bkmer, db_graph->kmer_size, bkstr);
  strbuf_sprintf(sbuf, "%s %zu\n", bkstr, subset->list.len);

  char orchar[2] = {0};
  orchar[FORWARD] = 'F';
  orchar[REVERSE] = 'R';
  const uint8_t *nseenptr;
  size_t klen;

  for(i = 0; i < subset->list.len; i++)
  {
    gpath = subset->list.data[i];
    nseenptr = gpath_set_get_nseen(gpset, gpath);
    klen = gpath_set_get_klen(gpset, gpath);
    strbuf_sprintf(sbuf, "%c %zu %u %u", orchar[gpath->orient], klen,
                                         gpath->num_juncs, (uint32_t)nseenptr[0]);

    for(col = 1; col < ncols; col++)
      strbuf_sprintf(sbuf, ",%u", (uint32_t)nseenptr[col]);

    strbuf_append_char(sbuf, ' ');
    strbuf_ensure_capacity(sbuf, sbuf->end + gpath->num_juncs + 2);
    binary_seq_to_str(gpath->seq, gpath->num_juncs, sbuf->b+sbuf->end);
    sbuf->end += gpath->num_juncs;
    strbuf_append_char(sbuf, '\n');
  }

  if(sbuf->end > DEFAULT_IO_BUFSIZE)
    _gpath_save_flush(gzout, sbuf, outlock);

  return 0; // => keep iterating
}

// Save paths for a single kmer
// @subset is temporary memory
void gpath_fwrite_single_kmer(hkey_t hkey, FILE *fout,
                              GPathSubset *subset,
                              const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;
  const size_t ncols = gpstore->gpset.ncols;
  GPath *first_gpath = gpath_store_fetch(gpstore, hkey);
  const GPath *gpath;
  char orchar[2] = {0};
  orchar[FORWARD] = 'F';
  orchar[REVERSE] = 'R';
  const uint8_t *nseenptr;
  size_t i, col, klen;

  // Load and sort paths for given kmer
  gpath_subset_reset(subset);
  gpath_subset_load_llist(subset, first_gpath);

  for(i = 0; i < subset->list.len; i++)
  {
    gpath = subset->list.data[i];
    nseenptr = gpath_set_get_nseen(gpset, gpath);
    klen = gpath_set_get_klen(gpset, gpath);
    fprintf(fout, "%c %zu %u %u", orchar[gpath->orient], klen,
                                  gpath->num_juncs, (uint32_t)nseenptr[0]);

    for(col = 1; col < ncols; col++)
      fprintf(fout, ",%u", (uint32_t)nseenptr[col]);

    fputc(' ', fout);
    binary_seq_print(gpath->seq, gpath->num_juncs, fout);
    fputc('\n', fout);
  }
}

typedef struct
{
  size_t threadid, nthreads;
  gzFile gzout;
  pthread_mutex_t *outlock;
  dBGraph *db_graph;
} GPathSaver;

static void gpath_save_thread(void *arg)
{
  GPathSaver wrkr = *(GPathSaver*)arg;
  const dBGraph *db_graph = wrkr.db_graph;

  GPathSubset subset;
  StrBuf sbuf;

  gpath_subset_alloc(&subset);
  gpath_subset_init(&subset, &wrkr.db_graph->gpstore.gpset);
  strbuf_alloc(&sbuf, 2 * DEFAULT_IO_BUFSIZE);

  HASH_ITERATE_PART(&db_graph->ht, wrkr.threadid, wrkr.nthreads,
                    _gpath_save_node,
                    wrkr.gzout, wrkr.outlock, &sbuf, &subset, db_graph);

  _gpath_save_flush(wrkr.gzout, &sbuf, wrkr.outlock);

  gpath_subset_dealloc(&subset);
  strbuf_dealloc(&sbuf);
}

// @hdrs is array of JSON headers of input files
void gpath_save(gzFile gzout, const char *path, size_t nthreads,
                cJSON **hdrs, size_t nhdrs,
                size_t *contig_len_hist, size_t hist_len,
                dBGraph *db_graph)
{
  ctx_assert(nthreads > 0);
  ctx_assert(gpath_set_has_nseen(&db_graph->gpstore.gpset));

  char npaths_str[50];
  ulong_to_str(db_graph->gpstore.num_paths, npaths_str);

  status("Saving %s paths to: %s", npaths_str, path);
  status("  using %zu threads", nthreads);

  // Write header
  _gpath_save_hdr(gzout, path, hdrs, nhdrs, contig_len_hist, hist_len, db_graph);

  // Print comments about the format
  gzputs(gzout, "\n");
  gzputs(gzout, "# This file was generated with Cortex\n");
  gzputs(gzout, "#   written by Isaac Turner <turner.isaac@gmail.com>\n");
  gzputs(gzout, "#   url: "CORTEX_URL"\n");
  gzputs(gzout, "# \n");
  gzputs(gzout, "# Comment lines begin with a # and are ignored, but must come after the header\n");
  gzputs(gzout, "# Format is:\n");
  gzputs(gzout, "#   [kmer] [num_paths] ...(ignored)\n");
  gzputs(gzout, "#   [FR] [num_kmers] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] ...(ignored)\n");
  gzputs(gzout, "\n");

  // Multithreaded
  GPathSaver *wrkrs = ctx_calloc(nthreads, sizeof(GPathSaver));
  pthread_mutex_t outlock;
  size_t i;

  if(pthread_mutex_init(&outlock, NULL) != 0) die("Mutex init failed");

  for(i = 0; i < nthreads; i++) {
    wrkrs[i] = (GPathSaver){.threadid = i,
                            .nthreads = nthreads,
                            .gzout = gzout,
                            .outlock = &outlock,
                            .db_graph = db_graph};
  }

  // Iterate over kmers writing paths
  util_run_threads(wrkrs, nthreads, sizeof(*wrkrs), nthreads, gpath_save_thread);

  pthread_mutex_destroy(&outlock);
  ctx_free(wrkrs);

  status("[GPathSave] Graph paths saved to %s", path);
}
