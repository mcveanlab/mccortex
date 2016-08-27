#include "global.h"
#include "commands.h"
#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "graphs_load.h"
#include "binary_kmer.h"

// TODO: add .ctp.gz indexing

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
  graph_file_open2(&gfile, ctx_path, "r+", true, 0);

  if(!file_filter_from_direct(&gfile.fltr))
    die("Cannot open graph file with a filter ('in.ctx:blah' syntax)");

  // Open output file
  FILE *fout = out_path ? futil_fopen_create(out_path, "w") : stdout;

  // Start
  size_t filencols = gfile.hdr.num_of_cols;
  size_t kmer_size = gfile.hdr.kmer_size;
  const char *path = file_filter_path(&gfile.fltr);

  size_t ncols = file_filter_into_ncols(&gfile.fltr);
  size_t kmer_mem = sizeof(BinaryKmer) + (sizeof(Edges)+sizeof(Covg))*filencols;

  if(block_size) {
    block_kmers = block_size / kmer_mem;
  } else if(!block_size && !block_kmers) {
    block_size = 4 * ONE_MEGABYTE;
    block_kmers = block_size / kmer_mem;
  }

  // Update block-size
  block_size = block_kmers * kmer_mem;

  status("[index] block bytes: %zu kmers: %zu; kmer bytes: %zu, hdr: %zu",
         block_size, block_kmers, kmer_mem, (size_t)gfile.hdr_size);

  if(block_kmers == 0) die("Cannot set block_kmers to zero");

  // Print header
  fputs("#block_start\tnext_block\tfirst_kmer\tkmer_idx\tnext_kmer_idx\n", fout);

  BinaryKmer bkmer = BINARY_KMER_ZERO_MACRO;
  BinaryKmer prev_bkmer = BINARY_KMER_ZERO_MACRO;
  Covg *covgs = ctx_malloc(ncols * sizeof(Covg));
  Edges *edges = ctx_malloc(ncols * sizeof(Edges));
  char bkmerstr[MAX_KMER_SIZE+1];

  size_t rem_block = block_size - kmer_mem; // block after first kmer
  char *tmp_mem = ctx_malloc(rem_block);

  // Read in file, print index
  size_t nblocks = 0;
  size_t bl_bytes = 0, bl_kmers = 0;
  size_t bl_byte_offset = gfile.hdr_size, bl_kmer_offset = 0;

  while(1)
  {
    if(!graph_file_read(&gfile, &bkmer, covgs, edges)) {
      status("Read kmer failed"); break; }
    binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
    if(nblocks > 0 && !binary_kmer_less_than(prev_bkmer,bkmer))
      die("File is not sorted: %s [%s]", bkmerstr, path);
    // We've already read one kmer entry, read rest of block
    bl_bytes = kmer_mem + gfr_fread_bytes(&gfile, tmp_mem, rem_block);
    bl_kmers = 1 + bl_bytes / kmer_mem;
    fprintf(fout, "%zu\t%zu\t%s\t%zu\t%zu\n",
            bl_byte_offset, bl_byte_offset+bl_bytes, bkmerstr,
            bl_kmer_offset, bl_kmer_offset+bl_kmers);
    bl_byte_offset += bl_bytes;
    bl_kmer_offset += bl_kmers;
    nblocks++;
    if(bl_kmers < block_kmers) {
      status("last block %zu < %zu; %zu vs %zu",
             bl_kmers, block_kmers, bl_bytes, block_size);
      break;
    }
    prev_bkmer = bkmer;
  }

  ctx_free(covgs);
  ctx_free(edges);
  ctx_free(tmp_mem);

  // done
  char num_kmers_str[50], num_blocks_str[50];
  char block_mem_str[50], block_kmers_str[50];
  ulong_to_str(bl_kmer_offset, num_kmers_str);
  ulong_to_str(nblocks, num_blocks_str);
  bytes_to_str(block_size, 1, block_mem_str);
  ulong_to_str(block_kmers, block_kmers_str);

  status("Read %s kmers in %s block%s (block size %s / %s kmers)",
         num_kmers_str, num_blocks_str, util_plural_str(nblocks),
         block_mem_str, block_kmers_str);

  if(fout != stdout) status("Saved to %s", out_path);

  graph_file_close(&gfile);
  fclose(fout);

  return EXIT_SUCCESS;
}
