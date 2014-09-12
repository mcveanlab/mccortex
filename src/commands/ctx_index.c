#include "global.h"
#include "commands.h"
#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "graph_format.h"
#include "binary_kmer.h"

// DEV: add .ctp.gz sorting

const char index_usage[] =
"usage: "CMD" index [options] <in.ctx>\n"
"\n"
"  Index a sorted cortex graph file (sort with `"CMD" sort` first).\n"
"\n"
"  -h, --help               This help message\n"
"  -q, --quiet              Silence status output normally printed to STDERR\n"
"  -f, --force              Overwrite output files\n"
"  -o, --out <out.ctx.idx>  Output file [default: STDOUT]\n"
"  -s, --block-size <S>     Block of <S> bytes [default: 4MB]\n"
"  -b, --block-kmers <B>    Block of <B> kmers\n"
"\n";

static struct option longopts[] =
{
  {"help",         no_argument,       NULL, 'h'},
  {"force",        no_argument,       NULL, 'f'},
  {"out",          required_argument, NULL, 'o'},
  {"block-size",   required_argument, NULL, 's'},
  {"block-kmers",  required_argument, NULL, 'b'},
  {NULL, 0, NULL, 0}
};

static inline size_t read_kmer(FILE *fin, const char *path,
                               size_t kmer_mem, char *tmp)
{
  size_t n = fread(tmp, 1, kmer_mem, fin);
  if(n == 0) return 0;
  if(n != kmer_mem) die("Invalid graph file: %s", path);
  return 1;
}

// Return number of kmers in the block
static inline size_t index_block(GraphFileReader *gfile,
                                 char *tmp, size_t kmer_mem,
                                 size_t block_kmers, size_t *offset,
                                 FILE *fout)
{
  size_t n_read, kmer_size = gfile->hdr.kmer_size;
  const char *path = file_filter_path(&gfile->fltr);

  BinaryKmer bkmer_start, bkmer_end;
  char kmer_start[MAX_KMER_SIZE+1], kmer_end[MAX_KMER_SIZE+1];

  if(!read_kmer(gfile->fh, path, kmer_mem, tmp)) return 0;

  memcpy(bkmer_start.b, tmp, sizeof(BinaryKmer));
  binary_kmer_to_str(bkmer_start, kmer_size, kmer_start);

  for(n_read = 1; n_read < block_kmers; n_read++) {
    if(!read_kmer(gfile->fh, path, kmer_mem, tmp))
      break;
  }

  memcpy(bkmer_end.b, tmp, sizeof(BinaryKmer));
  binary_kmer_to_str(bkmer_end, kmer_size, kmer_end);

  if(strcmp(kmer_start, kmer_end) >= 0)
    die("File is not sorted: %s vs %s [%s]", kmer_start, kmer_end, path);

  size_t block_start = *offset, block_size = n_read * kmer_mem;
  *offset += block_size;

  // Print index entry
  fprintf(fout, "%s %s %zu %zu %zu\n", kmer_start, kmer_end, n_read,
                                       block_start, block_size);

  return n_read;
}

int ctx_index(int argc, char **argv)
{
  const char *out_path = NULL;
  size_t block_size = 0, block_kmers = 0;

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
      case 'o': cmd_check(!out_path, cmd); out_path = optarg; break;
      case 'b':
        cmd_check(!block_kmers, cmd);
        block_kmers = cmd_size_nonzero(cmd, optarg);
        break;
      case 's':
        cmd_check(!block_size, cmd);
        block_size = cmd_size_nonzero(cmd, optarg);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" index -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(optind+1 != argc)
    cmd_print_usage("Require exactly one input graph file (.ctx)");

  if(block_size && block_kmers)
    cmd_print_usage("Cannot use --block-kmers and --block-size together");

  const char *ctx_path = argv[optind];

  //
  // Open Graph file
  //
  GraphFileReader gfile;
  memset(&gfile, 0, sizeof(GraphFileReader));
  graph_file_open2(&gfile, ctx_path, "r+");

  if(!file_filter_is_direct(&gfile.fltr))
    die("Cannot open graph file with a filter ('in.ctx:blah' syntax)");

  // Open output file
  FILE *fout = out_path ? futil_open_create(out_path, "w") : stdout;

  // Start
  size_t ncols = gfile.hdr.num_of_cols;
  size_t kmer_mem = sizeof(BinaryKmer) + (sizeof(Edges)+sizeof(Covg))*ncols;
  size_t nkmers, num_of_blocks = 0, num_kmers = 0, offset = 0;

  if(block_size) {
    block_kmers = block_size / kmer_mem;
  } else if(block_kmers) {
    block_size = block_kmers * kmer_mem;
  } else {
    block_size = 4 * ONE_MEGABYTE;
    block_kmers = block_size / kmer_mem;
  }

  if(block_kmers == 0) die("Cannot set block_kmers to zero");

  char tmp[kmer_mem];
  memset(tmp, 0, sizeof(tmp));

  // Read in file, print index
  while(1)
  {
    nkmers = index_block(&gfile, tmp, kmer_mem, block_kmers, &offset, fout);
    num_kmers += nkmers;
    num_of_blocks += (nkmers > 0);
    if(nkmers < block_kmers) break;
  }

  // done
  status("Read %zu kmers in %zu block%s", num_kmers, num_of_blocks,
                                          util_plural_str(num_of_blocks));

  if(fout != stdout) status("Saved to %s", out_path);

  graph_file_close(&gfile);
  fclose(fout);

  return EXIT_SUCCESS;
}
