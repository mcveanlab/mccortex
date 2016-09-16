#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
#include "graph_search.h"
#include "json_hdr.h"

const char server_usage[] =
"usage: "CMD" server [options] <in.ctx> [in2.ctx ...]\n"
"\n"
"  Interactively query the graph. Responds to STDOUT with JSON.\n"
"  Commands are:\n"
"  * 'info'     - print graph header\n"
"  * 'random'   - print a random kmer\n"
"  * 'ACACCAA'  - print information for the given kmer\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -p, --paths <in.ctp>  Load link file (can specify multiple times)\n"
"  -S, --single-line     Reponses on a single line\n"
"  -C, --coverages       Load coverages for kmers+links\n"
"  -E, --edges           Load per sample edges\n"
"  -D, --disk            Read from disk (one graph only, must be sorted)\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"paths",        required_argument, NULL, 'p'},
  {"single-line",  no_argument,       NULL, 'S'},
  {"coverages",    no_argument,       NULL, 'C'},
  {"edges",        no_argument,       NULL, 'E'},
  {"disk",         no_argument,       NULL, 'D'},
  {NULL, 0, NULL, 0}
};

#define MAX_RANDOM_TRIES 100

typedef struct {
  dBNode node;
  BinaryKmer bkey;
  Covg *covgs;
  Edges *edges;
  size_t ncols, nedges;
  bool binary_covgs;
} ServerQuery;

static void query_alloc(ServerQuery *q, size_t ncols,
                        bool binary_covgs, bool flatten_edges)
{
  q->covgs = ctx_calloc(ncols, sizeof(Covg));
  q->edges = ctx_calloc(ncols, sizeof(Edges));
  q->ncols = ncols;
  q->nedges = flatten_edges ? 1 : ncols;
  q->binary_covgs = binary_covgs;
}

static void query_dealloc(ServerQuery *q) {
  ctx_free(q->covgs);
  ctx_free(q->edges);
}

static inline void kmer_response(StrBuf *resp, ServerQuery q, bool pretty,
                                 const dBGraph *db_graph)
{
  size_t i;
  char keystr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(q.bkey, db_graph->kmer_size, keystr);

  strbuf_append_str(resp, "{");
  strbuf_append_str(resp, pretty ? "\n  " : " ");
  strbuf_append_str(resp, "\"key\": \"");
  strbuf_append_str(resp, keystr);
  strbuf_append_str(resp, "\", \"colours\": [");
  strbuf_append_ulong(resp, q.covgs[0]);
  for(i = 1; i < q.ncols; i++) {
    strbuf_append_char(resp, ',');
    strbuf_append_ulong(resp, q.covgs[i]);
  }
  strbuf_append_str(resp, "],");
  strbuf_append_str(resp, pretty ? "\n  " : " ");

  // Edges
  Edges uedges = 0; // get union of edges
  for(i = 0; i < q.nedges; i++) uedges |= q.edges[i];
  char edgesstr[9], left[5] = {0}, right[5] = {0}, *l = left, *r = right;
  db_node_get_edges_str(uedges, edgesstr);
  for(i = 0; i < 4; i++)
    if(edgesstr[i] != '.') { *l = toupper(edgesstr[i]); *(++l) = '\0'; }
  for(i = 4; i < 8; i++)
    if(edgesstr[i] != '.') { *r = edgesstr[i]; *(++r) = '\0'; }

  strbuf_append_str(resp, "\"left\": \"");
  strbuf_append_str(resp, left);
  strbuf_append_str(resp, "\", \"right\": \"");
  strbuf_append_str(resp, right);
  strbuf_append_str(resp, "\",");
  strbuf_append_str(resp, pretty ? "\n  " : " ");
  strbuf_append_str(resp, "\"edges\": \"");

  // Sample edges
  char sedges[3];
  for(i = 0; i < q.nedges; i++)
    strbuf_append_str(resp, edges_to_char(q.edges[i], sedges));

  strbuf_append_str(resp, "\",");
  strbuf_append_str(resp, pretty ? "\n  " : " ");
  strbuf_append_str(resp, "\"links\": [");

  // Links
  // {"forward": true, "juncs": "ACAA", "colours": [0,0,1]}
  size_t nlinks;
  const GPath *gpath = gpath_store_safe_fetch(&db_graph->gpstore, q.node.key);
  const GPathSet *gpset = &db_graph->gpstore.gpset;
  for(nlinks = 0; gpath != NULL; gpath = gpath->next, nlinks++)
  {
    if(nlinks) strbuf_append_str(resp, pretty ? ",\n            " : ", ");
    strbuf_append_str(resp, "{\"forward\": ");
    strbuf_append_str(resp, gpath->orient == FORWARD ? "true" : "false");
    strbuf_append_str(resp, ", \"juncs\": \"");

    // Print link sequence
    for(i = 0; i < gpath->num_juncs; i++)
      strbuf_append_char(resp, dna_nuc_to_char(binary_seq_get(gpath->seq, i)));

    // Print link colours
    // counts may be null if user did not specify -C,--coverages
    uint8_t *counts = gpath_set_get_nseen(gpset, gpath);
    strbuf_append_str(resp, "\", \"colours\": [");
    for(i = 0; i < db_graph->num_of_cols; i++) {
      if(i) strbuf_append_char(resp, ',');
      size_t count = counts ? counts[i]
                            : gpath_has_colour(gpath, gpset->ncols, i);
      strbuf_append_ulong(resp, count);
    }
    strbuf_append_str(resp, "]}");
  }

  strbuf_append_str(resp, pretty ? "]\n}\n" : "] }\n");
}

static inline void query_fetch_from_graph(ServerQuery *q, const dBGraph *db_graph)
{
  size_t i;
  for(i = 0; i < q->ncols; i++)
    q->covgs[i] = db_graph->col_covgs ? db_node_get_covg(db_graph, q->node.key, i)
                                      : db_node_has_col(db_graph, q->node.key, i);
  for(i = 0; i < q->nedges; i++)
    q->edges[i] = db_node_get_edges(db_graph, q->node.key, i);
}

static inline void query_fetch_from_disk(ServerQuery *q)
{
  size_t i;
  // Convert coverage to binary if required
  if(q->binary_covgs)
    for(i = 0; i < q->ncols; i++)
      q->covgs[i] = (q->covgs[i] > 0);
  // Flatten edges if we aren't outputting per sample edges
  if(q->nedges == 1)
    for(i = 1; i < q->ncols; i++)
      q->edges[0] |= q->edges[i];
}

/*
// Query: ACACCAA
{
  "key": "ACACAAA",
  "colours": [1,0,1],
  "left": "AC",
  "right": "",
  "edges": "c0",
  "links": [{"forward": true, "juncs": "ACAA", "colours": [0,0,1]},
            {"forward": false, "juncs": "TG", "colours": [1,0,0]}]
}
// Query: "ACCCCAC" (Not in graph)
{}
*/
/**
 * @param qstr    query string - must be "random" or kmer
 * @param resp    string buffer reset, then used to store response
 * @param pretty  pretty print JSON or one line JSON
 * @returns       true iff query was valid kmer
 */
static inline bool query_response(const char *qstr, ServerQuery q,
                                  StrBuf *resp, bool pretty,
                                  GraphFileSearch *disk, const dBGraph *db_graph)
{
  size_t qlen;
  strbuf_reset(resp);

  // query must be a kmer
  for(qlen = 0; qstr[qlen]; qlen++) {
    if(!char_is_acgt(qstr[qlen])) {
      strbuf_set(resp, "{\"error\": \"Invalid base\"}\n");
      return false;
    }
  }

  // Don't do anything if empty line
  if(qlen == 0) { return false; }

  if(qlen != db_graph->kmer_size) {
    strbuf_set(resp, "{\"error\": \"Doesn't match kmer size: ");
    strbuf_append_ulong(resp, db_graph->kmer_size);
    strbuf_append_str(resp, "\"}\n");
    return false;
  }

  BinaryKmer bkmer = binary_kmer_from_str(qstr, db_graph->kmer_size);
  q.bkey = binary_kmer_get_key(bkmer, db_graph->kmer_size);
  q.node.orient = (binary_kmer_eq(bkmer, q.bkey) ? FORWARD : REVERSE);

  if(disk == NULL) {
    // Fetch from graph
    q.node.key = hash_table_find(&db_graph->ht, q.bkey);
    if(q.node.key == HASH_NOT_FOUND) { strbuf_set(resp, "{}\n"); return true; }
    query_fetch_from_graph(&q, db_graph);
  }
  else {
    if(!graph_search_find(disk, q.bkey, q.covgs, q.edges)) {
      strbuf_set(resp, "{}\n");
      return true;
    }
    query_fetch_from_disk(&q);
  }

  kmer_response(resp, q, pretty, db_graph);
  return true;
}

// Reply with a random kmer
static inline void request_random(ServerQuery q, StrBuf *resp, bool pretty,
                                  GraphFileSearch *disk, const dBGraph *db_graph)
{
  strbuf_reset(resp);
  if(disk == NULL) {
    q.node.key = db_graph_rand_node(db_graph, MAX_RANDOM_TRIES);
    if(q.node.key == HASH_NOT_FOUND) { strbuf_set(resp, "{}\n"); return; }
    q.node.orient = FORWARD;
    q.bkey = db_node_get_bkey(db_graph, q.node.key);
    query_fetch_from_graph(&q, db_graph);
  }
  else {
    graph_search_rand(disk, &q.bkey, q.covgs, q.edges);
    query_fetch_from_disk(&q);
  }
  kmer_response(resp, q, pretty, db_graph);
}

static char* make_info_json_str(cJSON **hdrs, size_t nhdrs,
                                bool pretty, size_t nkmers,
                                const dBGraph *db_graph)
{
  cJSON *json = cJSON_CreateObject();
  json_hdr_make_std(json, NULL, hdrs, nhdrs, db_graph, nkmers);

  cJSON *paths = cJSON_CreateObject();
  cJSON_AddItemToObject(json, "paths", paths);

  // Add command specific header fields
  const GPathStore *gpstore = &db_graph->gpstore;
  cJSON_AddNumberToObject(paths, "num_kmers_with_paths", gpstore->num_kmers_with_paths);
  cJSON_AddNumberToObject(paths, "num_paths", gpstore->num_paths);
  cJSON_AddNumberToObject(paths, "path_bytes", gpstore->path_bytes);

  char *info_txt = pretty ? cJSON_Print(json) : cJSON_PrintUnformatted(json);
  cJSON_Delete(json);

  return info_txt;
}

int ctx_server(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

  bool pretty = true;
  bool binary_covgs = true; // Binary coverage instead of full coverage
  bool per_col_edges = false; // Load per sample or pooled edges
  bool use_disk = false;

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_push(&gpfiles, &tmp_gpfile, 1);
        break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'S': cmd_check(pretty, cmd); pretty = false; break;
      case 'C': cmd_check(binary_covgs, cmd); binary_covgs = false; break;
      case 'E': cmd_check(!per_col_edges, cmd); per_col_edges = true; break;
      case 'D': cmd_check(!use_disk, cmd); use_disk = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" server -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(optind >= argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  const size_t num_gfiles = argc - optind;
  char **graph_paths = argv + optind;

  ctx_assert(num_gfiles > 0);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t i, ncols;
  size_t ctx_max_kmers = 0, ctx_sum_kmers = 0;
  size_t ctp_max_kmers = 0, ctp_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  gpath_reader_count_kmers(gpfiles.b, gpfiles.len, &ctp_max_kmers, &ctp_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, gpfiles.b, gpfiles.len, -1);

  if(use_disk && num_gfiles > 1)
    cmd_print_usage("Can only use --disk with one sorted graph file");

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem = 0, path_mem = 0;

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +

  if(use_disk && gpfiles.len == 0)
  {
    kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                          memargs.mem_to_use_set,
                                          memargs.num_kmers,
                                          memargs.num_kmers_set,
                                          0, 0, 0, false, &graph_mem);
  }
  else
  {
    bits_per_kmer = sizeof(BinaryKmer)*8 + // kmer
                    sizeof(Edges)*8 * (per_col_edges ? ncols : 1) + // edges
                    (binary_covgs ? 1 : sizeof(Covg)*8) * ncols + // covgs
                    (gpfiles.len > 0 ? sizeof(GPath*)*8 : 0); // links

    kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                          memargs.mem_to_use_set,
                                          memargs.num_kmers,
                                          memargs.num_kmers_set,
                                          bits_per_kmer,
                                          use_disk ? ctp_max_kmers : ctx_max_kmers,
                                          use_disk ? ctp_sum_kmers : ctx_sum_kmers,
                                          false, &graph_mem);

    if(gpfiles.len)
    {
      // Paths memory
      size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
      path_mem = gpath_reader_mem_req(gpfiles.b, gpfiles.len,
                                      ncols, rem_mem,
                                      !binary_covgs); // load path counts

      // Shift path store memory from graphs->paths
      graph_mem -= sizeof(GPath*)*kmers_in_hash;
      path_mem  += sizeof(GPath*)*kmers_in_hash;
      cmd_print_mem(path_mem, "paths");
    }
  }

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  // Allocate memory
  int allocflags = DBG_ALLOC_EDGES | (binary_covgs ? DBG_ALLOC_NODE_IN_COL
                                                   : DBG_ALLOC_COVGS);
  if(use_disk) allocflags = 0;

  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size,
                 ncols, per_col_edges ? ncols : 1, kmers_in_hash,
                 allocflags);

  // Paths - allocates nothing if gpfiles.len == 0
  gpath_reader_alloc_gpstore(gpfiles.b, gpfiles.len,
                             path_mem, !binary_covgs,
                             &db_graph);

  //
  // Load graphs
  //
  GraphFileSearch *disk = NULL;

  if(use_disk) {
    // Only load graph info
    graph_load_ginfo(&db_graph, &gfiles[0]);
    disk = graph_search_new(&gfiles[0]);
  }
  else {
    GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);
    gprefs.empty_colours = true;
    for(i = 0; i < num_gfiles; i++) {
      graph_load(&gfiles[i], gprefs, NULL);
      graph_file_close(&gfiles[i]);
      gprefs.empty_colours = false;
    }
  }

  // Load link files
  int link_flags = use_disk ? GPATH_ADD_MISSING_KMERS : GPATH_DIE_MISSING_KMERS;
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_load(&gpfiles.b[i], link_flags, &db_graph);

  hash_table_print_stats(&db_graph.ht);

  // Create array of cJSON** from input files
  cJSON **hdrs = ctx_malloc(gpfiles.len * sizeof(cJSON*));
  for(i = 0; i < gpfiles.len; i++) hdrs[i] = gpfiles.b[i].json;

  // Construct cJSON
  size_t nkmers_in_graph = use_disk ? ctx_max_kmers : hash_table_nkmers(&db_graph.ht);
  char *info_txt = make_info_json_str(hdrs, gpfiles.len, pretty,
                                      nkmers_in_graph, &db_graph);
  ctx_free(hdrs);

  // Close input link files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_close(&gpfiles.b[i]);
  gpfile_buf_dealloc(&gpfiles);

  // Answer queries
  StrBuf line, response;
  strbuf_alloc(&line, 1024);
  strbuf_alloc(&response, 1024);
  size_t nqueries = 0, nbad_queries = 0;
  bool success;

  ServerQuery q;
  query_alloc(&q, db_graph.num_of_cols, binary_covgs, !per_col_edges);

  // Read from input
  while(1)
  {
    fprintf(stdout, "> "); fflush(stdout);
    if(futil_fcheck(strbuf_reset_readline(&line, stdin), stdin, "STDIN") == 0) {
      fprintf(stdout, "\n");
      break;
    }
    strbuf_chomp(&line);
    if(strcasecmp(line.b,"q") == 0 || strcasecmp(line.b,"quit") == 0) { break; }
    else if(strcasecmp(line.b,"info") == 0) {
      fputs(info_txt, stdout);
      fputc('\n', stdout);
      fflush(stdout);
    }
    else if(strcasecmp(line.b,"random") == 0) {
      request_random(q, &response, pretty, disk, &db_graph);
      fputs(response.b, stdout);
      fflush(stdout);
    }
    else {
      success = query_response(line.b, q, &response, pretty, disk, &db_graph);
      if(response.end) {
        fputs(response.b, stdout);
        fflush(stdout);
      }
      nbad_queries += !success;
    }
    nqueries += (line.end > 0);
  }

  char nstr[50], badstr[50];
  ulong_to_str(nqueries, nstr);
  ulong_to_str(nbad_queries, badstr);
  status("Answered %s queries, %s bad queries", nstr, badstr);

  query_dealloc(&q);

  if(disk) {
    graph_search_destroy(disk);
    graph_file_close(&gfiles[0]);
  }
  ctx_free(gfiles);

  free(info_txt);
  strbuf_dealloc(&line);
  strbuf_dealloc(&response);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
