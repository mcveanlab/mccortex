#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "gpath_reader.h"
#include "gpath_checks.h"

const char server_usage[] =
"usage: "CMD" server [options] <in.ctx> [in2.ctx ...]\n"
"\n"
"  Interactively query the graph. Responds to STDOUT with JSON.\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -p, --paths <in.ctp>  Load path file (can specify multiple times)\n"
// "  -P,--port 80          Run as HTTP server on given port\n"
// "  -C,--coverages        Print per sample coverages"
// "  -E,--edges            Print per sample edges"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"paths",        required_argument, NULL, 'p'},
  {NULL, 0, NULL, 0}
};

/*
// Query: ACACCAA
{
  "key": "ACACAAA",
  "colours": [1,0,1],
  "left": "AC",
  "right": "",
  "links": [{"forward": true, "juncs": "ACAA", "colours": [0,0,1]},
            {"forward": false, "juncs": "TG", "colours": [1,0,0]}]
}
// Query: ACCCCAC
{}
*/
// @returns true iff query was valid
static inline bool print_kmer_json(const char *kstr, const dBGraph *db_graph)
{
  FILE *fout = stdout;
  size_t klen, i, col;
  dBNode node;
  char keystr[MAX_KMER_SIZE+1], *ptr;

  for(klen = 0; kstr[klen]; klen++) {
    if(!char_is_acgt(kstr[klen])) {
      fprintf(fout, "{\"error\": \"Invalid base\"}\n");
      return false;
    }
  }

  if(klen == 0) { return false; }
  if(klen != db_graph->kmer_size) {
    fprintf(fout, "{\"error\": \"Doesn't match kmer size: %zu\"}\n",
            db_graph->kmer_size);
    return false;
  }

  node = db_graph_find_str(db_graph, kstr);
  if(node.key == HASH_NOT_FOUND) {
    fprintf(fout, "{}\n");
    return true;
  }

  // Get upper case kmer key
  memcpy(keystr, kstr, klen+1);
  for(ptr = keystr; *ptr; ptr++) *ptr = toupper(*ptr);
  if(node.orient == REVERSE) dna_reverse_complement_str(keystr, klen);

  fprintf(fout, "{\n"
                "  \"key\": \"%s\", \"colours\": [%c",
          keystr, '0' + db_node_has_col(db_graph, node.key, 0));

  // 0/1 in colours
  char colstr[3] = ",0";
  for(col = 1; col < db_graph->num_of_cols; col++) {
    colstr[1] = '0'+db_node_has_col(db_graph, node.key, 0);
    fputs(colstr, fout);
  }

  // Edges
  Edges edges = db_node_get_edges_union(db_graph, node.key);    
  char edgesstr[9], left[5], right[5], *l, *r;
  db_node_get_edges_str(edges, edgesstr);
  for(l = left, i = 0; i < 4; i++)
    if(edgesstr[i] != '.') { *l = toupper(edgesstr[i]); *(++l) = '\0'; }
  for(r = right, i = 4; i < 8; i++)
    if(edgesstr[i] != '.') { *r = edgesstr[i]; *(++r) = '\0'; }

  fprintf(fout, "], \"left\": \"%s\", \"right\": \"%s\",\n"
                "  \"links\": [", left, right);

  // Links
  // {"forward": true, "juncs": "ACAA", "colours": [0,0,1]}
  size_t nlinks;
  GPath *gpath = gpath_store_fetch(&db_graph->gpstore, node.key);
  for(nlinks = 0; gpath != NULL; gpath = gpath->next, nlinks++)
  {
    if(nlinks) fputs(",\n            ", fout);
    fprintf(fout, "{\"forward\": %s, \"juncs\": \"",
                  gpath->orient == FORWARD ? "true" : "false");

    // Print link sequence
    for(i = 0; i < gpath->num_juncs; i++)
      fputc(dna_nuc_to_char(binary_seq_get(gpath->seq, i)), fout);

    // Print link colours
    fprintf(fout, "\", \"colours\": [%c",
            '0'+gpath_has_colour(gpath, db_graph->num_of_cols, 0));
    for(col = 1; col < db_graph->num_of_cols; col++) {
      colstr[1] = '0'+gpath_has_colour(gpath, db_graph->num_of_cols, col);
      fputs(colstr, fout);
    }
    fputs("]}", fout);
  }

  fputs("]\n}\n", fout);

  return true;
}

int ctx_server(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;

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
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem;

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +

  bits_per_kmer = sizeof(BinaryKmer)*8 + sizeof(Edges)*8 +
                  (gpfiles.len > 0 ? sizeof(GPath*)*8 : 0) +
                  ncols;

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        ctx_max_kmers, ctx_sum_kmers,
                                        false, &graph_mem);

  // Paths memory
  size_t rem_mem = memargs.mem_to_use - MIN2(memargs.mem_to_use, graph_mem);
  path_mem = gpath_reader_mem_req(gpfiles.b, gpfiles.len, ncols, rem_mem, false);

  // Shift path store memory from graphs->paths
  graph_mem -= sizeof(GPath*)*kmers_in_hash;
  path_mem  += sizeof(GPath*)*kmers_in_hash;
  cmd_print_mem(path_mem, "paths");

  size_t total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash,
                 DBG_ALLOC_EDGES | DBG_ALLOC_NODE_IN_COL);

  // Paths
  gpath_reader_alloc_gpstore(gpfiles.b, gpfiles.len, path_mem, false, &db_graph);

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

  // Close input path files
  for(i = 0; i < gpfiles.len; i++)
    gpath_reader_close(&gpfiles.b[i]);
  gpfile_buf_dealloc(&gpfiles);

  // Answer queries
  StrBuf line;
  strbuf_alloc(&line, 1024);
  size_t nqueries = 0;

  // Read from input
  while(strbuf_reset_readline(&line, stdin))
  {
    strbuf_chomp(&line);
    if(print_kmer_json(line.b, &db_graph))
      nqueries++;
  }

  status("Answered %zu queries", nqueries);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}


// GCGTTACAATATCGTATTGGGTTCGTGACCAACACTCCCATTTCTTGATATGACGCCATCAAACGACTACACGGAGACCGGCCGAGCATGGCAACCCGCACGACTGCATCATCTCCATCAATCCAACCATACTCCCGGACTTACCCCTGCCCCGGGCGCAGCAGTCCTTAAGATCAGGAACTGGGGTGTACGACGGCCTCGCTGACACGGTACCAGCCGTGCACCGATGCTGCTAGGCACCCGTCGCCTGCTCAAGAAATGGCTGGGTTCAATAAGCGTTTGTGAGTGCTTCGACTCGTTAGGATGTAATTAGGGCCAGTAGTCAACCAGCGCTAGTGAGAATATGATAGAGATTTCGCAAAGTCCTTGGTATACAGGATCTCAACCCACAGACTGCGGAGGCTGTGGTGCCATCATCGGACTCACTACGTCCTTGTCAGGCCTAACCTTTCAGGGCGGCAAGCTACGGTTACCTGACCGAAGTCTTATTCACAGTTCGGTAGCTCCAATCATTGCGAGGTTAGCTTAACGCCTGACATTACCTGGCAAACATGCTCCTTTCACGACCGTTTATCGGCGCGATTTGATATCCACTTG
