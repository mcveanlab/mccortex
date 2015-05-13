#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "gpath_reader.h"
#include "gpath_checks.h"
#include "gpath_save.h"
#include "json_hdr.h"

const char pview_usage[] =
"usage: "CMD" pview [options] [-p <in.ctp>] <in.ctx> [in2.ctx ...]\n"
"\n"
"  View cortex path files (.ctp).\n"
"\n"
"  -h, --help             This help message\n"
"  -q, --quiet            Silence status output normally printed to STDERR\n"
// "  -f, --force            Overwrite output files\n"
// "  -o, --out <out.txt>    Output file [required]\n"
"  -m, --memory <mem>     Memory to use\n"
"  -n, --nkmers <kmers>   Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -p, --paths <in.ctp>   Load path file (can specify multiple times)\n"
// "  -H, --header-only      Only print the header (no paths)\n"
// "  -P, --paths-only       Only print the paths (no header)\n"
"\n";


static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"paths",        required_argument, NULL, 'p'},
  {"force",        no_argument,       NULL, 'f'},
// command specific
  {"header-only",  no_argument,       NULL, 'H'},
  {"paths-only",   no_argument,       NULL, 'P'},
  {NULL, 0, NULL, 0}
};

static int _print_paths(hkey_t hkey,
                        StrBuf *sbuf, GPathSubset *subset,
                        dBNodeBuffer *nbuf, SizeBuffer *jposbuf,
                        FILE *fout, const dBGraph *db_graph)
{
  strbuf_reset(sbuf);
  gpath_save_sbuf(hkey, sbuf, subset, nbuf, jposbuf, db_graph);
  fputs(sbuf->b, fout);
  return 0; // 0 => keep iterating
}

static cJSON* _get_header(GPathFileBuffer *gpfiles, const dBGraph *db_graph)
{
  // Load contig hist distribution
  size_t i;
  ZeroSizeBuffer *contig_histgrms = ctx_calloc(db_graph->num_of_cols,
                                               sizeof(ZeroSizeBuffer));

  for(i = 0; i < db_graph->num_of_cols; i++)
    zsize_buf_alloc(&contig_histgrms[i], 512);

  size_t j, fromcol, intocol;
  for(i = 0; i < gpfiles->len; i++) {
    for(j = 0; j < file_filter_num(&gpfiles->b[i].fltr); j++) {
      fromcol = file_filter_fromcol(&gpfiles->b[i].fltr, j);
      intocol = file_filter_intocol(&gpfiles->b[i].fltr, j);
      gpath_reader_load_contig_hist(gpfiles->b[i].json,
                                    gpfiles->b[i].fltr.path.b,
                                    fromcol, &contig_histgrms[intocol]);
    }
  }

  cJSON *hdrs[gpfiles->len];
  for(i = 0; i < gpfiles->len; i++) hdrs[i] = gpfiles->b[i].json;
  cJSON *json = gpath_save_mkhdr("STDOUT", NULL, NULL, hdrs, gpfiles->len,
                                 contig_histgrms, db_graph->num_of_cols,
                                 db_graph);

  for(i = 0; i < db_graph->num_of_cols; i++)
    zsize_buf_dealloc(&contig_histgrms[i]);

  ctx_free(contig_histgrms);

  return json;
}

int ctx_pview(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  bool header_only = false, paths_only = false;

  GPathReader tmp_gpfile;
  GPathFileBuffer gpfiles;
  gpfile_buf_alloc(&gpfiles, 8);

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg);
        gpfile_buf_push(&gpfiles, &tmp_gpfile, 1);
        break;
      case 'H': cmd_check(!header_only, cmd); header_only = true; break;
      case 'P': cmd_check(!paths_only,  cmd); paths_only  = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" clean -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";

  if(optind >= argc)   cmd_print_usage("Please give input graph files");
  if(gpfiles.len == 0) cmd_print_usage("Please give input path files");

  if(header_only && paths_only) cmd_print_usage("Cannot use both -H and -P");

  // Use remaining args as graph files
  char **gfile_paths = argv + optind;
  size_t i, num_gfiles = (size_t)(argc - optind);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(gfile_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  // Check graph + paths are compatible
  graphs_gpaths_compatible(gfiles, num_gfiles, gpfiles.b, gpfiles.len, -1);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 +
                  sizeof(GPath*)*8 + ncols;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.b, gpfiles.len, ncols, rem_mem, true);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Open output file
  //
  FILE *fout = futil_fopen_create(out_path, "w");

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.b, gpfiles.len, path_mem, true, &db_graph);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Load path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_load(&gpfiles.b[i], GPATH_DIE_MISSING_KMERS, &db_graph);

  // Generate merged header
  if(!paths_only) {
    cJSON *json = _get_header(&gpfiles, &db_graph);
    json_hdr_fprint(json, fout);
    fputs(ctp_explanation_comment, fout);
    cJSON_Delete(json);
  }

  if(!header_only)
  {
    // Print paths
    StrBuf sbuf;
    dBNodeBuffer nbuf;
    SizeBuffer jposbuf;
    GPathSubset subset;

    strbuf_alloc(&sbuf, 4096);
    db_node_buf_alloc(&nbuf, 1024);
    size_buf_alloc(&jposbuf, 256);
    gpath_subset_alloc(&subset);
    gpath_subset_init(&subset, &db_graph.gpstore.gpset);

    HASH_ITERATE(&db_graph.ht, _print_paths,
                 &sbuf, &subset, &nbuf, &jposbuf,
                 fout, &db_graph);

    strbuf_dealloc(&sbuf);
    db_node_buf_dealloc(&nbuf);
    size_buf_dealloc(&jposbuf);
    gpath_subset_dealloc(&subset);
  }

  if(fout != stdout) fclose(fout);

  // Close input path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_close(&gpfiles.b[i]);
  gpfile_buf_dealloc(&gpfiles);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
