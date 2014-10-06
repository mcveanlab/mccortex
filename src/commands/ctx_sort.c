#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "graph_format.h"
#include "binary_kmer.h"

// DEV: add .ctp.gz sorting

const char sort_usage[] =
"usage: "CMD" sort [options] <in.ctx>\n"
"\n"
"  Sort a cortex graph file.\n"
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

static inline int cmp_bkmer(const void *aa, const void *bb)
{
  BinaryKmer b1, b2;
  const char *a = *(const char *const*)aa, *b = *(const char *const*)bb;
  memcpy(b1.b, a, sizeof(BinaryKmer));
  memcpy(b2.b, b, sizeof(BinaryKmer));
  return binary_kmers_cmp(b1, b2);
}

// Sort either graph file or paths file. Pointers must point to binary kmer or
// test representation of kmer
static inline void sort_block(char **entries, size_t num)
{
  qsort(entries, num, sizeof(char*), cmp_bkmer);
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
  graph_file_open2(&gfile, ctx_path, out_path ? "r" : "r+", 0);

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
  FILE *fout = out_path ? futil_open_create(out_path, "w") : NULL;

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

  // Read in file
  // This fseek not needed
  // if(fseek(gfile.fh, gfile.hdr_size, SEEK_SET) == -1) die("fseek failed");
  size_t nkread = fread(mem, kmer_mem, num_kmers, gfile.fh);

  // check we are at the end of the file
  int b;
  if(nkread == num_kmers && (b = fgetc(gfile.fh)) != -1) {
    die("More kmers in file than believed (%i; kmers: %zu ncols: %zu).",
        b, num_kmers, ncols);
  }

  num_kmers = nkread;

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
    if(fseek(gfile.fh, gfile.hdr_size, SEEK_SET) == -1) die("fseek failed");
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
