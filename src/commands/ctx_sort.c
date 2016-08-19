#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "graphs_load.h"
#include "graph_writer.h"
#include "binary_kmer.h"

// TODO: add .ctp.gz sorting

const char sort_usage[] =
"usage: "CMD" sort [options] <in.ctx>\n"
"\n"
"  Sort a cortex graph file. Loads entire graph into memory then sorts.\n"
"\n"
"  -h, --help              This help message\n"
"  -q, --quiet             Silence status output normally printed to STDERR\n"
"  -f, --force             Overwrite output files\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -o, --out <out.ctx>     Output file [default: overwrite input]\n"
"\n";

static struct option longopts[] =
{
  {"help",         no_argument,       NULL, 'h'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"out",          required_argument, NULL, 'o'},
  {NULL, 0, NULL, 0}
};

// Sort graph file entries. Pointers must point to binary kmer
static inline void sort_block(char **entries, size_t num)
{
  qsort(entries, num, sizeof(char*), binary_kmers_qcmp_unaligned_ptrs);
}

int ctx_sort(int argc, char **argv)
{
  const char *out_path = NULL;
  struct MemArgs memargs = MEM_ARGS_INIT;

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
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" sort -h` for help. Bad option: %s", argv[optind-1]);
      default: die("Bad option: [%c]: %s", c, cmd);
    }
  }

  if(optind+1 != argc)
    cmd_print_usage("Require exactly one input graph file (.ctx)");

  const char *ctx_path = argv[optind];

  //
  // Open Graph file
  //
  GraphFileReader gfile;
  memset(&gfile, 0, sizeof(GraphFileReader));
  graph_file_open2(&gfile, ctx_path, out_path ? "r" : "r+", true, 0);

  if(!file_filter_is_direct(&gfile.fltr))
    die("Cannot open graph file with a filter ('in.ctx:blah' syntax)");

  size_t num_kmers, memory;

  // Reading from a stream
  if(gfile.num_of_kmers < 0) {
    if(!memargs.num_kmers_set)
      die("If reading from a stream, must give -n <num_kmers>");
    num_kmers = memargs.num_kmers;
  }
  else num_kmers = gfile.num_of_kmers;

  // Open output path (if given)
  FILE *fout = out_path ? futil_fopen_create(out_path, "w") : NULL;

  size_t i;
  size_t ncols = gfile.hdr.num_of_cols;
  size_t kmer_mem = sizeof(BinaryKmer) + (sizeof(Edges)+sizeof(Covg))*ncols;

  memory = (sizeof(char*) + kmer_mem) * num_kmers;

  char mem_str[50];
  bytes_to_str(memory, 1, mem_str);

  if(memory > memargs.mem_to_use)
    die("Require at least %s memory", mem_str);

  status("[memory] Total: %s", mem_str);

  char *mem = ctx_malloc(kmer_mem * num_kmers);
  char **kmers = ctx_malloc(num_kmers*sizeof(char*));

  // Read in whole file
  // if(graph_file_fseek(gfile, gfile.hdr_size, SEEK_SET) != 0) die("fseek failed");
  size_t nkread = gfr_fread_bytes(&gfile, mem, num_kmers*kmer_mem);

  if(nkread != num_kmers*kmer_mem)
    die("Could only read %zu bytes [<%zu]", nkread, num_kmers*kmer_mem);

  // check we are at the end of the file
  char tmpc;
  if(gfr_fread_bytes(&gfile, &tmpc, 1) != 0) {
    die("More kmers in file than believed (kmers: %zu ncols: %zu).",
        num_kmers, ncols);
  }

  status("Read %zu kmers with %zu colour%s", num_kmers,
         ncols, util_plural_str(ncols));

  for(i = 0; i < num_kmers; i++)
    kmers[i] = mem + kmer_mem*i;

  sort_block(kmers, num_kmers);

  // Print
  if(out_path != NULL) {
    // saving to a different destination - write header
    graph_write_header(fout, &gfile.hdr);
  }
  else {
    // Directly manipulating gfile.fh here, using it to write later
    // Not doing any more reading
    if(fseek(gfile.fh, gfile.hdr_size, SEEK_SET) != 0) die("fseek failed");
    fout = gfile.fh;
  }

  for(i = 0; i < num_kmers; i++)
    if(fwrite(kmers[i], 1, kmer_mem, fout) != kmer_mem)
      die("Cannot write to file");

  if(out_path) fclose(fout);

  graph_file_close(&gfile);
  ctx_free(kmers);
  ctx_free(mem);

  return EXIT_SUCCESS;
}
