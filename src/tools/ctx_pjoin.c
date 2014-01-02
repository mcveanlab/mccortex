#include "global.h"
#include "string_buffer.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_store.h"
#include "path_format.h"

static const char usage[] =
"usage: "CMD" pjoin [options] <out.ctp> <in0.ctp> [offset:]in1.ctp[:0,2-4] ...\n"
"  Merge cortex path files.\n"
"\n"
" Options:\n"
"   -m <mem>       Memory to use (required) recommend 80G for human\n"
"   -n <kmers>     Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"   -f <in.ctx>    Get number of hash table entries from graph file\n"
"   --overlap      Merge corresponding colours from each graph file\n"
"   --flatten      Dump into a single colour graph\n"
"   --outcols <C>  How many 'colours' should the output file have\n"
"\n"
"  Files can be specified with specific colours: samples.ctp:2,3\n"
"  Offset specifies where to load the first colour: 3:samples.ctp\n";

int ctx_pjoin(CmdArgs *args)
{
  cmd_accept_options(args, "mnf", usage);
  // cmd_require_options(args, "m", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  // if(!args->num_kmers_set && !args->input_file_set)
  //   print_usage(usage, "Please specify -n <num-kmers> or -f <in.ctx>");
  // if(args->num_kmers_set && args->input_file_set)
  //   print_usage(usage, "Please specify only ONE of -n <num-kmers> or -f <in.ctx>");

  boolean overlap = false, flatten = false;
  size_t output_ncols = 0;
  int argi;

  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++) {
    if(strcasecmp(argv[argi],"--overlap") == 0) {
      if(overlap) warn("overlap specified twice");
      overlap = true;
    }
    else if(strcasecmp(argv[argi],"--flatten") == 0) {
      if(flatten) warn("flatten specified twice");
      flatten = true;
    }
    else if(strcasecmp(argv[argi],"--outcols") == 0) {
      if(argi+1 == argc || !parse_entire_size(argv[argi+1], &output_ncols) ||
         output_ncols == 0) {
        print_usage(usage, "--outcols <C> needs an integer argument > 0");
      }
      argi++;
    }
    else {
      print_usage(usage, "Unknown argument '%s'", argv[argi]);
    }
  }

  if(argc - argi < 2)
    print_usage(usage, "Please specify output and input paths");

  const char *out_ctp_path = argv[argi++];

  // argi .. argend-1 are graphs to load
  size_t num_pfiles = argc - argi;
  char **paths = argv + argi;

  //
  // Open all path files
  //
  size_t i, ncols, max_cols = 0, sum_cols = 0, total_cols;
  size_t ctp_max_path_kmers = 0, ctp_max_path_bytes = 0;
  PathFileReader pfiles[num_pfiles];

  for(i = 0; i < num_pfiles; i++)
  {
    pfiles[i] = INIT_PATH_READER;
    path_file_open(&pfiles[i], paths[i], true);

    if(pfiles[0].hdr.kmer_size != pfiles[i].hdr.kmer_size) {
      print_usage(usage, "Kmer sizes don't match [%u vs %u]",
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
    ctp_max_path_bytes = MAX2(ctp_max_path_bytes, pfiles[i].hdr.num_path_bytes);

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
    print_usage(usage, "You specified --outcols %zu but inputs need at %zu colours",
                output_ncols, total_cols);
  }

  // Open graph file to get number of kmers is passed
  uint64_t num_kmers = ctp_max_path_kmers;
  GraphFileReader gfile = INIT_GRAPH_READER;

  if(args->num_kmers_set) {
    num_kmers = MAX2(num_kmers, args->num_kmers);
  }

  if(args->input_file_set) {
    graph_file_open(&gfile, args->input_file, true);
    if(gfile.hdr.kmer_size != pfiles[0].hdr.kmer_size) {
      warn("Kmer-sizes don't match graph: %u paths: %u [graph: %s path: %s]",
           gfile.hdr.kmer_size, pfiles[0].hdr.kmer_size,
           gfile.fltr.file_path.buff, pfiles[0].fltr.file_path.buff);
    }
    num_kmers = MAX2(num_kmers, gfile.hdr.num_of_kmers);
  }

  if(args->num_kmers_set && args->num_kmers < num_kmers) {
    char num_kmers_str[100], args_num_kmers_str[100];
    ulong_to_str(num_kmers, num_kmers_str);
    ulong_to_str(args->num_kmers, args_num_kmers_str);
    warn("Using %s kmers instead of (-n) %s", num_kmers_str, args_num_kmers_str);
  }

  // if(num_kmers < ctp_max_path_kmers) {
  //   print_usage(usage, "Please set a larger -n <kmers> (needs to be > %zu)",
  //               ctp_max_path_kmers);
  // }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, total_mem;
  char path_mem_str[100];

  // Each kmer stores a pointer to its list of paths
  bits_per_kmer = sizeof(uint64_t)*8;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, num_kmers,
                                        false, &graph_mem);

  // Path Memory
  size_t tmppathsize, req_path_mem;
  tmppathsize = paths_merge_needs_tmp(pfiles, num_pfiles) ? ctp_max_path_bytes : 0;
  req_path_mem = tmppathsize + ctp_max_path_bytes;

  size_t req_mem = graph_mem + req_path_mem;
  char req_mem_str[100];
  bytes_to_str(req_mem, 1, req_mem_str);

  path_mem = args->mem_to_use - graph_mem;
  bytes_to_str(path_mem, 1, path_mem_str);
  status("[memory] paths: %s\n", path_mem_str);

  if(path_mem < req_path_mem)
    die("Not enough memory - require %s. Decrease -n or increase -m", req_mem_str);

  total_mem = graph_mem + path_mem;
  cmd_check_mem_limit(args, total_mem);

  // Set up graph and PathStore
  dBGraph db_graph;

  db_graph_alloc(&db_graph, pfiles[0].hdr.kmer_size, output_ncols, 0, kmers_in_hash);

  for(i = 0; i < num_pfiles; i++)
    path_file_set_graph_sample_names(&pfiles[i], &db_graph);

  db_graph.kmer_paths = malloc2(db_graph.ht.capacity * sizeof(uint64_t));
  memset((void*)db_graph.kmer_paths, 0xff, db_graph.ht.capacity * sizeof(uint64_t));

  uint8_t *path_store = malloc2(path_mem);
  path_store_init(&db_graph.pdata, path_store, path_mem, output_ncols);

  // Temorary memory to load paths into
  uint8_t *tmppdata = tmppathsize > 0 ? malloc2(tmppathsize) : NULL;

  // Open output file
  FILE *fout = fopen(out_ctp_path, "w");
  if(fout == NULL) die("Cannot open output file: %s", out_ctp_path);

  //
  // Set up file header
  //
  PathFileHeader pheader = INIT_PATH_FILE_HDR;
  pheader.version = CTX_PATH_FILEFORMAT;
  pheader.kmer_size = pfiles[0].hdr.kmer_size;
  pheader.num_of_cols = output_ncols;

  paths_header_alloc(&pheader, output_ncols);

  // Load path files
  boolean add_kmers = true;
  paths_format_merge(pfiles, num_pfiles, add_kmers, tmppdata, tmppathsize, &db_graph);

  for(i = 0; i < num_pfiles; i++)
    path_file_set_header_sample_names(&pfiles[i], &pheader);

  // Dump paths file
  setvbuf(fout, NULL, _IOFBF, CTP_BUF_SIZE);
  paths_header_update(&pheader, &db_graph.pdata);
  paths_format_write_header(&pheader, fout);
  paths_format_write_optimised_paths(&db_graph, fout);
  fclose(fout);

  char pnum_str[100], pbytes_str[100], pkmers_str[100];
  ulong_to_str(pheader.num_of_paths, pnum_str);
  bytes_to_str(pheader.num_path_bytes, 1, pbytes_str);
  ulong_to_str(pheader.num_kmers_with_paths, pkmers_str);

  status("Paths written to: %s\n", out_ctp_path);
  status("  %s paths, %s path-bytes, %s kmers", pnum_str, pbytes_str, pkmers_str);

  free((void *)db_graph.kmer_paths);
  free(path_store);
  if(tmppdata != NULL) free(tmppdata);

  db_graph_dealloc(&db_graph);

  graph_file_dealloc(&gfile);
  paths_header_dealloc(&pheader);
  for(i = 0; i < num_pfiles; i++) path_file_dealloc(&pfiles[i]);

  return EXIT_SUCCESS;
}
