#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_format.h"
#include "binary_kmer.h"
#include "supernode.h"

static const char usage[] =
"usage: "CMD" supernodes [options] <out.fa.gz> <in.ctx> [<in2.ctx> ...]\n"
"  Print supernodes with k-1 bases of overlap.\n"
// "  Prints to stdout, messages to stderr\n"
"  Options:  -m <mem> | -h <kmers> | --fastg\n";

static size_t ncap, supernode_idx = 0;
static hkey_t *nodes;
static Orientation *orients;

// In fastg, supernodes are edges between vertices
// we compute vertices by finding the lowest kmer key at a junction and
// converting it to hexidecimal
static inline void get_vertex(dBGraph *db_graph, hkey_t node, Orientation orient,
                              char *str, boolean first)
{
  BinaryKmer bkmer = db_graph_oriented_bkmer(db_graph, node, orient);
  if(first) binary_kmer_set_last_nuc(&bkmer, 0);
  else binary_kmer_left_shift_one_base(&bkmer, db_graph->kmer_size);
  bkmer = db_node_get_key(bkmer, db_graph->kmer_size-1);
  binary_kmer_to_hex(bkmer, db_graph->kmer_size-1, str);
}

static void dump_supernodes(hkey_t node, gzFile gzout, dBGraph *db_graph,
                            uint64_t *visited, boolean use_fastg)
{
  boolean cycle;
  size_t i, len;
  char fastg_from[MAX_KMER_SIZE+1], fastg_to[MAX_KMER_SIZE+1];

  if(!bitset_has(visited, node))
  {
    len = supernode_find(db_graph, node, &nodes, &orients, &cycle, &ncap);
    for(i = 0; i < len; i++) bitset_set(visited, nodes[i]);
    if(use_fastg) {
      get_vertex(db_graph, nodes[0], !orients[0], fastg_from, true);
      get_vertex(db_graph, nodes[len-1], orients[len-1], fastg_to, false);
      gzprintf(gzout, ">%s:%s;\n", fastg_from, fastg_to);
    }
    else gzprintf(gzout, ">supernode%zu\n", supernode_idx);
    supernode_gzprint(gzout, db_graph, nodes, orients, len);
    gzputc(gzout, '\n');
    supernode_idx++;
  }
}

// Returns 0 on success, otherwise != 0
int ctx_supernodes(CmdArgs *args)
{
  cmd_accept_options(args, "mp");

  int argi, argc = args->argc;
  char **argv = args->argv;
  if(argc < 2) print_usage(usage, NULL);

  boolean use_fastg = false;
  size_t i, num_files;
  char **paths;
  uint64_t max_ctx_kmers = 0;

  if(strcasecmp(argv[0],"--fastg") == 0) {
    use_fastg = true;
    argv++;
    argc--;
  }

  for(argi = 0; argi < argc; argi++)
    if(argv[argi][0] == '-')
      print_usage(usage, "Unexpected option: %s", argv[argi]);

  if(argc < 2) print_usage(usage, NULL);

  const char *out_path = argv[0];
  argv++;
  argc--;

  if(!test_file_writable(out_path))
    print_usage(usage, "Cannot write to output file: %s", out_path);

  num_files = argc;
  paths = argv;
  GraphFileReader files[num_files];

  for(i = 0; i < num_files; i++)
  {
    files[i] = INIT_GRAPH_READER;
    int ret = graph_file_open(&files[i], paths[i], false);

    if(ret == 0)
      print_usage(usage, "Cannot read input graph file: %s", paths[i]);
    else if(ret < 0)
      print_usage(usage, "Input graph file isn't valid: %s", paths[i]);

    if(files[0].hdr.kmer_size != files[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
                  files[0].hdr.kmer_size, files[i].hdr.kmer_size);
    }

    max_ctx_kmers = MAX2(max_ctx_kmers, files[i].hdr.num_of_kmers);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;
  bits_per_kmer = sizeof(Edges) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, max_ctx_kmers, true);

  status("Output in %s format\n", use_fastg ? "FASTG" : "FASTA");

  //
  // Open output file
  //
  gzFile gzout = gzopen(out_path, "w");
  if(gzout == NULL) die("Cannot open output file: %s", out_path);

  //
  // Allocate memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, files[0].hdr.kmer_size, 1, 1, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity, sizeof(Edges));

  // Visited
  size_t numwords64 = round_bits_to_words64(db_graph.ht.capacity);
  uint64_t *visited = calloc2(numwords64, sizeof(uint64_t));

  SeqLoadingPrefs prefs = {.db_graph = &db_graph,
                           .boolean_covgs = true,
                           .must_exist_in_graph = false,
                           .empty_colours = false};

  for(i = 0; i < num_files; i++) {
    files[i].flatten = true;
    files[i].intocol = 0;
    graph_load(&files[i], &prefs, NULL);
  }

  if(use_fastg)
  {
    gzprintf(gzout, "#FASTG:begin;\n#FASTG:version=1.0:assembly_name=\"%s\";\n",
                    db_graph.ginfo[0].sample_name.buff);
  }

  hash_table_print_stats(&db_graph.ht);

  ncap = 2048;
  nodes = malloc(ncap * sizeof(hkey_t));
  orients = malloc(ncap * sizeof(Orientation));

  // dump supernodes
  status("Writing supernodes to %s\n", out_path);
  HASH_TRAVERSE(&db_graph.ht, dump_supernodes, gzout, &db_graph, visited, use_fastg);

  status("Dumped %zu supernodes\n", supernode_idx);

  gzclose(gzout);

  free(nodes);
  free(orients);
  free(visited);
  free(db_graph.col_edges);
  db_graph_dealloc(&db_graph);

  for(i = 0; i < num_files; i++) graph_file_dealloc(&files[i]);

  return EXIT_SUCCESS;
}
