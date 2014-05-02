#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_store.h"
#include "path_format.h"

const char pjoin_usage[] =
"usage: "CMD" pjoin [options] <in1.ctp> [[offset:]in2.ctp[:0,2-4] ...]\n"
"  Merge cortex path files.\n"
"\n"
" Options:\n"
"   -o <out.ctp>     Output file [required]\n"
"   -m <mem>         Memory to use (required) recommend 80G for human\n"
"   -n <kmers>       Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   --graph <in.ctx> Get number of hash table entries from graph file\n"
"   --overlap        Merge corresponding colours from each graph file\n"
"   --flatten        Dump into a single colour graph\n"
"   --outcols <C>    How many 'colours' should the output file have\n"
"   --noredundant    Remove redundant paths\n"
"\n"
"  Files can be specified with specific colours: samples.ctp:2,3\n"
"  Offset specifies where to load the first colour: 3:samples.ctp\n";

int ctx_pjoin(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  // Have already checked we have at least 1 argument

  int argi;
  bool overlap = false, flatten = false, noredundant = false;
  size_t output_ncols = 0;
  char *graph_file = NULL;

  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++) {
    if(strcmp(argv[argi],"--graph") == 0) {
      if(argi+1 == argc) cmd_print_usage("--graph <in.ctx> needs an argument");
      graph_file = argv[argi+1];
      argi++;
    }
    else if(strcmp(argv[argi],"--overlap") == 0) {
      if(overlap) warn("overlap specified twice");
      overlap = true;
    }
    else if(strcmp(argv[argi],"--flatten") == 0) {
      if(flatten) warn("flatten specified twice");
      flatten = true;
    }
    else if(strcmp(argv[argi],"--noredundant") == 0) {
      if(noredundant) warn("noredundant specified twice");
      noredundant = true;
    }
    else if(strcmp(argv[argi],"--outcols") == 0) {
      if(argi+1 == argc || !parse_entire_size(argv[argi+1], &output_ncols) ||
         output_ncols == 0) {
        cmd_print_usage("--outcols <C> needs an integer argument > 0");
      }
      argi++;
    }
    else {
      cmd_print_usage("Unknown argument '%s'", argv[argi]);
    }
  }

  if(argi >= argc)
    cmd_print_usage("Please specify at least one input file");

  const char *out_ctp_path = args->output_file;

  // argi .. argend-1 are graphs to load
  size_t num_pfiles = (size_t)(argc - argi);
  char **paths = argv + argi;

  //
  // Open all path files
  //
  size_t i, ncols, max_cols = 0, sum_cols = 0, total_cols;
  size_t ctp_max_path_kmers = 0;
  PathFileReader pfiles[num_pfiles];

  for(i = 0; i < num_pfiles; i++)
  {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], paths[i], true);

    if(pfiles[0].hdr.kmer_size != pfiles[i].hdr.kmer_size) {
      cmd_print_usage("Kmer sizes don't match [%u vs %u]",
                  pfiles[0].hdr.kmer_size, pfiles[i].hdr.kmer_size);
    }

    if(flatten) {
      pfiles[i].fltr.flatten = true;
      file_filter_update_intocol(&pfiles[i].fltr, 0);
    }

    ncols = path_file_usedcols(&pfiles[i]);
    max_cols = MAX2(max_cols, ncols);
    sum_cols += ncols;
    ctp_max_path_kmers = MAX2(ctp_max_path_kmers, pfiles[i].hdr.num_kmers_with_paths);

    file_filter_status(&pfiles[i].fltr);
  }

  if(flatten) total_cols = 1;
  else if(overlap) total_cols = max_cols;
  else {
    total_cols = 0;
    for(i = 0; i < num_pfiles; i++) {
      size_t offset = total_cols;
      total_cols += path_file_usedcols(&pfiles[i]);
      file_filter_update_intocol(&pfiles[i].fltr, pfiles[i].fltr.intocol+offset);
    }
  }

  if(output_ncols == 0) output_ncols = total_cols;
  else if(total_cols > output_ncols) {
    cmd_print_usage("You specified --outcols %zu but inputs need at %zu colours",
                output_ncols, total_cols);
  }

  // Open graph file to get number of kmers is passed
  uint64_t num_kmers = ctp_max_path_kmers;
  GraphFileReader gfile = INIT_GRAPH_READER;

  if(args->num_kmers_set) {
    num_kmers = MAX2(num_kmers, args->num_kmers);
  }

  if(graph_file != NULL) {
    graph_file_open(&gfile, graph_file, true);
    if(gfile.hdr.kmer_size != pfiles[0].hdr.kmer_size) {
      warn("Kmer-sizes don't match graph: %u paths: %u [graph: %s path: %s]",
           gfile.hdr.kmer_size, pfiles[0].hdr.kmer_size,
           gfile.fltr.file_path.buff, pfiles[0].fltr.file_path.buff);
    }
    num_kmers = MAX2(num_kmers, gfile.num_of_kmers);
  }

  if(args->num_kmers_set && args->num_kmers < num_kmers) {
    char num_kmers_str[100], args_num_kmers_str[100];
    ulong_to_str(num_kmers, num_kmers_str);
    ulong_to_str(args->num_kmers, args_num_kmers_str);
    warn("Using %s kmers instead of (-n) %s", num_kmers_str, args_num_kmers_str);
  }

  // if(num_kmers < ctp_max_path_kmers) {
  //   cmd_print_usage("Please set a larger -n <kmers> (needs to be > %zu)",
  //               ctp_max_path_kmers);
  // }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;
  size_t path_mem_req, path_mem, total_mem;

  // Each kmer stores a pointer to its list of paths
  bits_per_kmer = sizeof(uint64_t)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        num_kmers, num_kmers,
                                        false, &graph_mem);

  // Path Memory
  path_mem_req = path_files_mem_required(pfiles, num_pfiles, false, false);
  path_mem = MAX2(args->mem_to_use - graph_mem, path_mem_req);
  cmd_print_mem(path_mem, "paths");

  total_mem = graph_mem + path_mem;

  cmd_check_mem_limit(args, total_mem);

  // Set up graph and PathStore
  dBGraph db_graph;
  db_graph_alloc(&db_graph, pfiles[0].hdr.kmer_size, output_ncols, 0, kmers_in_hash);

  for(i = 0; i < num_pfiles; i++)
    path_file_set_graph_sample_names(&pfiles[i], &db_graph);

  path_store_alloc(&db_graph.pstore, path_mem, false,
                   db_graph.ht.capacity, output_ncols);

  // Open output file
  FILE *fout = fopen(out_ctp_path, "w");
  if(fout == NULL) die("Cannot open output file: %s", out_ctp_path);

  //
  // Set up file header
  //
  PathFileHeader pheader = INIT_PATH_FILE_HDR;
  pheader.version = CTX_PATH_FILEFORMAT;
  pheader.kmer_size = pfiles[0].hdr.kmer_size;
  pheader.num_of_cols = (uint32_t)output_ncols;

  paths_header_alloc(&pheader, output_ncols);

  // Load path files
  bool add_kmers = true;
  paths_format_merge(pfiles, num_pfiles, add_kmers,
                     noredundant, args->max_work_threads, &db_graph);

  for(i = 0; i < num_pfiles; i++)
    path_file_set_header_sample_names(&pfiles[i], &pheader);

  // Dump paths file
  setvbuf(fout, NULL, _IOFBF, CTP_BUF_SIZE);
  paths_header_update(&pheader, &db_graph.pstore);
  paths_format_write_header(&pheader, fout);
  paths_format_write_optimised_paths(&db_graph, fout);
  fclose(fout);

  char pnum_str[100], pbytes_str[100], pkmers_str[100];
  ulong_to_str(pheader.num_of_paths, pnum_str);
  bytes_to_str(pheader.num_path_bytes, 1, pbytes_str);
  ulong_to_str(pheader.num_kmers_with_paths, pkmers_str);

  status("Paths written to: %s\n", out_ctp_path);
  status("  %s paths, %s path-bytes, %s kmers", pnum_str, pbytes_str, pkmers_str);

  path_store_dealloc(&db_graph.pstore);
  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&gfile);
  paths_header_dealloc(&pheader);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  return EXIT_SUCCESS;
}
