#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graphs_load.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
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
"  -p, --paths <in.ctp>  Load path file (can specify multiple times)\n"
"  -S, --single-line     Reponses on a single line\n"
"  -C, --coverages       Load per sample coverages\n"
"  -E, --edges           Load per sample edges\n"
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
  {NULL, 0, NULL, 0}
};

#define MAX_RANDOM_TRIES 100

static inline void kmer_response(StrBuf *resp, dBNode node, const char *keystr,
                                 bool pretty, const dBGraph *db_graph)
{
  size_t i, col;

  strbuf_append_str(resp, "{");
  strbuf_append_str(resp, pretty ? "\n  " : " ");
  strbuf_append_str(resp, "\"key\": \"");
  strbuf_append_str(resp, keystr);
  strbuf_append_str(resp, "\", \"colours\": [");
  for(col = 0; col < db_graph->num_of_cols; col++) {
    if(col) strbuf_append_char(resp, ',');
    Covg covg = db_graph->col_covgs ? db_node_get_covg(db_graph, node.key, col)
                                    : db_node_has_col(db_graph, node.key, col);
    strbuf_append_ulong(resp, covg);
  }
  strbuf_append_str(resp, "],");
  strbuf_append_str(resp, pretty ? "\n  " : " ");

  // Edges
  Edges edges = db_node_get_edges_union(db_graph, node.key);
  char edgesstr[9], left[5] = {0}, right[5] = {0}, *l = left, *r = right;
  db_node_get_edges_str(edges, edgesstr);
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
  for(i = 0; i < db_graph->num_edge_cols; i++) {
    edges_to_char(db_node_get_edges(db_graph, node.key, i), sedges);
    strbuf_append_str(resp, sedges);
  }

  strbuf_append_str(resp, "\",");
  strbuf_append_str(resp, pretty ? "\n  " : " ");
  strbuf_append_str(resp, "\"links\": [");

  // Links
  // {"forward": true, "juncs": "ACAA", "colours": [0,0,1]}
  size_t nlinks;
  const GPath *gpath = gpath_store_safe_fetch(&db_graph->gpstore, node.key);
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
    for(col = 0; col < db_graph->num_of_cols; col++) {
      if(col) strbuf_append_char(resp, ',');
      size_t count = counts ? counts[col]
                            : gpath_has_colour(gpath, gpset->ncols, col);
      strbuf_append_ulong(resp, count);
    }
    strbuf_append_str(resp, "]}");
  }

  strbuf_append_str(resp, pretty ? "]\n}\n" : "] }\n");
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
static inline bool query_response(const char *qstr, StrBuf *resp, bool pretty,
                                  const dBGraph *db_graph)
{
  size_t qlen;
  dBNode node;
  char keystr[MAX_KMER_SIZE+1], *ptr;
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

  node = db_graph_find_str(db_graph, qstr);
  if(node.key == HASH_NOT_FOUND) {
    strbuf_set(resp, "{}\n");
    return true;
  }

  // Get upper case kmer key
  memcpy(keystr, qstr, qlen+1);
  for(ptr = keystr; *ptr; ptr++) *ptr = toupper(*ptr);
  if(node.orient == REVERSE) dna_reverse_complement_str(keystr, qlen);
  kmer_response(resp, node, keystr, pretty, db_graph);
  return true;
}

// Reply with a random kmer
static inline void request_random(StrBuf *resp, bool pretty,
                                  const dBGraph *db_graph)
{
  dBNode node;
  char keystr[MAX_KMER_SIZE+1];
  strbuf_reset(resp);

  hkey_t hkey = db_graph_rand_node(db_graph, MAX_RANDOM_TRIES);
  if(hkey == HASH_NOT_FOUND) { strbuf_set(resp, "{}\n"); }
  node.key = hkey;
  node.orient = FORWARD;
  BinaryKmer bkmer = db_node_get_bkmer(db_graph, node.key);
  binary_kmer_to_str(bkmer, db_graph->kmer_size, keystr);
  kmer_response(resp, node, keystr, pretty, db_graph);
}

static char* make_info_json_str(cJSON **hdrs, size_t nhdrs,
                                bool pretty, const dBGraph *db_graph)
{
  cJSON *json = cJSON_CreateObject();
  json_hdr_make_std(json, NULL, hdrs, nhdrs, db_graph);

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
  // Per sample coverage and edges
  bool load_covgs = false, load_edges = false;

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
      case 'C': cmd_check(!load_covgs, cmd); load_covgs = true; break;
      case 'E': cmd_check(!load_edges, cmd); load_edges = true; break;
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
  size_t i, ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, gpfiles.b, gpfiles.len, -1);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem = 0;

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +

  bits_per_kmer = sizeof(BinaryKmer)*8 + // kmer
                  sizeof(Edges)*8 * (load_edges ? ncols : 1) + // edges
                  sizeof(Covg)*8 * (load_covgs ? ncols : 0) + // covgs
                  (gpfiles.len > 0 ? sizeof(GPath*)*8 : 0) + // links
                  ncols; // in colour

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, &graph_mem);

  if(gpfiles.len)
  {
    // Paths memory
    size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
    path_mem = gpath_reader_mem_req(gpfiles.b, gpfiles.len,
                                    ncols, rem_mem,
                                    load_covgs); // load path counts

    // Shift path store memory from graphs->paths
    graph_mem -= sizeof(GPath*)*kmers_in_hash;
    path_mem  += sizeof(GPath*)*kmers_in_hash;
    cmd_print_mem(path_mem, "paths");
  }

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size,
                 ncols, load_edges ? ncols : 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL |
                   (load_covgs ? DBG_ALLOC_COVGS : 0));

  // Paths - allocates nothing if gpfiles.len == 0
  gpath_reader_alloc_gpstore(gpfiles.b, gpfiles.len,
                             path_mem, load_covgs,
                             &db_graph);

  //
  // Load graphs
  //
  GraphLoadingPrefs gprefs = graph_loading_prefs(&db_graph);
  gprefs.empty_colours = true;

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, NULL);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_load(&gpfiles.b[i], GPATH_DIE_MISSING_KMERS, &db_graph);

  // Create array of cJSON** from input files
  cJSON **hdrs = ctx_malloc(gpfiles.len * sizeof(cJSON*));
  for(i = 0; i < gpfiles.len; i++) hdrs[i] = gpfiles.b[i].json;

  // Construct cJSON
  char *info_txt = make_info_json_str(hdrs, gpfiles.len, pretty, &db_graph);
  ctx_free(hdrs);

  // Close input path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_close(&gpfiles.b[i]);
  gpfile_buf_dealloc(&gpfiles);

  // Answer queries
  StrBuf line, response;
  strbuf_alloc(&line, 1024);
  strbuf_alloc(&response, 1024);
  size_t nqueries = 0, nbad_queries = 0;
  bool success;

  // Read from input
  while(1)
  {
    fprintf(stdout, "> "); fflush(stdout);
    if(futil_fcheck(strbuf_reset_readline(&line, stdin), stdin, "STDIN") == 0)
      break;
    strbuf_chomp(&line);
    if(strcmp(line.b,"q") == 0) { break; }
    else if(strcmp(line.b,"info") == 0) {
      fputs(info_txt, stdout);
      fputc('\n', stdout);
      fflush(stdout);
    }
    else if(strcmp(line.b,"random") == 0) {
      request_random(&response, pretty, &db_graph);
      fputs(response.b, stdout);
      fflush(stdout);
    }
    else {
      success = query_response(line.b, &response, pretty, &db_graph);
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

  free(info_txt);
  strbuf_dealloc(&line);
  strbuf_dealloc(&response);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
