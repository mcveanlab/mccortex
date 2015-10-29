#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "graphs_load.h"
#include "hash_mem.h" // for calculating mem usage

const char view_usage[] =
"usage: "CMD" view [options] <in.ctx>\n"
"\n"
"  View a cortex graph as a list of kmers with coverage and edges\n"
"\n"
"  -h, --help   This help message\n"
"  -q, --quiet  Silence status output normally printed to STDERR\n"
"\n"
"  -k, --kmers  Print kmers\n"
"  -c, --check  Check kmers\n"
"  -i, --info   Print info\n"
"\n"
" Default is [--info --check]\n"
"\n";

int print_info = 0, parse_kmers = 0, print_kmers = 0;

static struct option longopts[] =
{
  {"help",  no_argument, NULL,        'h'},
  {"kmers", no_argument, &print_kmers,  1},
  {"check", no_argument, &parse_kmers,  1},
  {"info",  no_argument, &print_info,   1},
  {NULL, 0, NULL, 0}
};

static void print_header(GraphFileHeader *h, size_t num_of_kmers)
{
  printf("version: %u\n", h->version);
  printf("kmer size: %u\n", h->kmer_size);
  printf("bitfields: %u\n", h->num_of_bitfields);
  printf("colours: %u\n", h->num_of_cols);

  char num_kmers_str[50];
  ulong_to_str(num_of_kmers, num_kmers_str);
  printf("number of kmers: %s\n", num_kmers_str);
  printf("----\n");

  size_t i;
  for(i = 0; i < h->num_of_cols; i++)
  {
    GraphInfo *ginfo = h->ginfo + i;

    printf("Colour %zu:\n", i);

    if(h->version >= 6)
    {
      // Version 6 only output
      printf("  sample name: '%s'\n", ginfo->sample_name.b);
    }

    char total_sequence_str[100];
    ulong_to_str(ginfo->total_sequence, total_sequence_str);

    printf("  mean input contig length: %u\n", ginfo->mean_read_length);
    printf("  total sequence loaded:    %s\n", total_sequence_str);

    if(h->version >= 6)
    {
      // Version 6 only output
      printf("  sequence error rate: %Lf\n", ginfo->seq_err);

      ErrorCleaning *ec = &ginfo->cleaning;
      printf("  tip clipping: %s\n", (ec->cleaned_tips == 0 ? "no" : "yes"));

      printf("  remove low coverage supernodes: %s [threshold: <%u]\n",
             ec->cleaned_snodes ? "yes" : "no",
             ec->clean_snodes_thresh);

      printf("  remove low coverage kmers: %s [threshold: <%u]\n",
             ec->cleaned_kmers ? "yes" : "no",
             ec->clean_kmers_thresh);

      printf("  cleaned against graph: %s [against: '%s']\n",
             ec->is_graph_intersection ? "yes" : "no",
             ec->intersection_name.b);
    }
  }
}

#define loading_warning(fmt,...) { num_warnings++; warn(fmt, ##__VA_ARGS__);}
#define loading_error(fmt,...) { num_errors++; warn(fmt, ##__VA_ARGS__);}

int ctx_view(int argc, char **argv)
{
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
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        cmd_print_usage("`"CMD" view -h` for help. Bad option: %s", argv[optind-1]);
      default:
        cmd_print_usage("Programmer fail. Tell Isaac.");
    }
  }

  if(print_kmers) parse_kmers = 1;

  if(!print_info && !parse_kmers && !print_kmers)
    print_info = parse_kmers = 1;

  if(optind+1 != argc) cmd_print_usage("Require one input graph file (.ctx)");

  char *path = argv[optind];
  size_t num_errors = 0, num_warnings = 0;

  GraphFileReader gfile;
  memset(&gfile, 0, sizeof(gfile));
  int ret = graph_file_open(&gfile, path);
  if(ret == 0) die("Cannot open file: %s", path);

  if(print_info)
  {
    char fsize_str[50];
    bytes_to_str((size_t)gfile.file_size, 0, fsize_str);
    printf("Loading file: %s\n", file_filter_path(&gfile.fltr));
    printf("File size: %s\n", fsize_str);
    printf("----\n");
  }

  size_t i, col, ncols = file_filter_into_ncols(&gfile.fltr);
  size_t kmer_size = gfile.hdr.kmer_size;
  ctx_assert(ncols > 0);

  GraphFileHeader hdr;
  memset(&hdr, 0, sizeof(hdr));
  graph_file_merge_header(&hdr, &gfile);

  uint64_t nkmers_read = 0, nkmers_loaded = 0;
  uint64_t num_all_zero_kmers = 0, num_zero_covg_kmers = 0;
  uint64_t *col_nkmers, *col_sum_covgs;
  col_nkmers = ctx_calloc(ncols, sizeof(col_nkmers[0]));
  col_sum_covgs = ctx_calloc(ncols, sizeof(col_sum_covgs[0]));

  // Print header
  if(print_info) print_header(&hdr, gfile.num_of_kmers);

  BinaryKmer bkmer;
  Covg covgs[ncols], keep_kmer;
  Edges edges[ncols];

  bool direct_read = file_filter_is_direct(&gfile.fltr);

  if(parse_kmers || print_kmers)
  {
    if(print_info && print_kmers) printf("----\n");

    for(; graph_file_read_reset(&gfile, &bkmer, covgs, edges); nkmers_read++)
    {
      // If kmer has no covg in any samples -> don't load
      keep_kmer = 0;
      for(col = 0; col < ncols; col++) {
        col_nkmers[col] += (covgs[col] > 0);
        col_sum_covgs[col] += covgs[col];
        keep_kmer |= covgs[col];
      }

      if(!direct_read && !keep_kmer) continue;
      nkmers_loaded++;

      /* Kmer Checks */
      // graph_file_read_reset() already checks for:
      // 1. oversized kmers
      // 2. kmers with covg 0 in all colours
      // 3. edges without coverage in a colour

      // Check for all-zeros (i.e. all As kmer: AAAAAA)
      uint64_t kmer_words_or = 0;

      for(i = 0; i < hdr.num_of_bitfields; i++)
        kmer_words_or |= bkmer.b[i];

      if(kmer_words_or == 0)
      {
        if(num_all_zero_kmers == 1)
        {
          loading_error("more than one all 'A's kmers seen [index: %"PRIu64"]\n",
                        nkmers_read);
        }

        num_all_zero_kmers++;
      }

      // Check covg is 0 for all colours
      for(i = 0; i < ncols && covgs[i] == 0; i++);
      num_zero_covg_kmers += (i == ncols);

      // Print
      if(print_kmers)
        db_graph_print_kmer2(bkmer, covgs, edges, ncols, kmer_size, stdout);
    }
  }

  // check for various reading errors
  if(errno != 0)
    loading_error("errno set [%i]: %s\n", (int)errno, strerror(errno));

  int err = ferror(gfile.fh);
  if(err != 0)
    loading_error("occurred after file reading [%i]\n", err);

  char nstr[50];

  if(print_kmers || parse_kmers)
  {
    // file_size is set to -1 if we are reading from a stream,
    // therefore won't be able to check number of kmers read
    if(gfile.file_size != -1 && nkmers_read != (uint64_t)gfile.num_of_kmers) {
      loading_warning("Expected %zu kmers, read %zu\n",
                      (size_t)gfile.num_of_kmers, (size_t)nkmers_read);
    }

    if(num_all_zero_kmers > 1)
    {
      loading_error("%s all-zero-kmers seen\n",
                    ulong_to_str(num_all_zero_kmers, nstr));
    }

    if(num_zero_covg_kmers > 0)
    {
      loading_warning("%s kmers have no coverage in any colour\n",
                      ulong_to_str(num_zero_covg_kmers, nstr));
    }
  }

  // Count warnings printed by graph_file_reader.c
  num_warnings += gfile.error_zero_covg;
  num_warnings += gfile.error_missing_covg;

  // Can only print these stats if we're read in the kmers
  if((print_kmers || parse_kmers) && print_info)
  {
    // print kmer coverage per sample
    printf("\n---- Per colour stats\n");
    printf("num. kmers:");
    for(col = 0; col < ncols; col++)
      printf("\t%s", ulong_to_str(col_nkmers[col], nstr));
    printf("\n");
    printf("sum coverage:");
    for(col = 0; col < ncols; col++)
      printf("\t%s", ulong_to_str(col_sum_covgs[col], nstr));
    printf("\n");
    printf("kmer coverage:");
    for(col = 0; col < ncols; col++)
      printf("\t%.2f", safe_frac(col_sum_covgs[col], col_nkmers[col]));
    printf("\n");

    // Overall stats
    uint64_t sum_covgs = 0;
    double mean_kmer_covg = 0.0;
    for(col = 0; col < ncols; col++) sum_covgs += col_sum_covgs[col];
    mean_kmer_covg = nkmers_loaded ? (double)sum_covgs / nkmers_loaded : 0.0;

    printf("\n---- Overall stats\n");
    printf("Total kmers:    %s\n", ulong_to_str(nkmers_loaded, nstr));
    printf("Total coverage: %s\n", ulong_to_str(sum_covgs, nstr));
    printf("Mean coverage:  %s\n", double_to_str(mean_kmer_covg, 2, nstr));
  }

  if(print_info)
  {
    // Print memory stats
    uint64_t mem, capacity, num_buckets, req_capacity;
    uint8_t bucket_size;

    req_capacity = (size_t)(gfile.num_of_kmers / IDEAL_OCCUPANCY);
    capacity = hash_table_cap(req_capacity, &num_buckets, &bucket_size);
    mem = ht_mem(bucket_size, num_buckets,
                 sizeof(BinaryKmer)*8 + ncols*(sizeof(Covg)+sizeof(Edges))*8);

    char memstr[100], capacitystr[100], bucket_size_str[100], num_buckets_str[100];
    bytes_to_str(mem, 1, memstr);
    ulong_to_str(capacity, capacitystr);
    ulong_to_str(bucket_size, bucket_size_str);
    ulong_to_str(num_buckets, num_buckets_str);

    size_t mem_height = (size_t)__builtin_ctzl(num_buckets);

    printf("\n---- Memory\n");
    printf("memory required: %s [capacity: %s]\n", memstr, capacitystr);
    printf("  bucket size: %s; number of buckets: %s\n",
            bucket_size_str, num_buckets_str);
    printf("  --kmer_size %zu --mem_height %zu --mem_width %i\n",
           kmer_size, mem_height, bucket_size);
  }

  if((print_kmers || parse_kmers) && print_info)
  {
    printf("\n----\n");
    if(num_warnings > 0 || num_errors > 0) {
      printf("Warnings: %zu; Errors: %zu\n",
              (size_t)num_warnings, (size_t)num_errors);
    }
    if(num_errors == 0)
      printf(num_warnings ? "Graph may be ok\n" : "Graph is valid\n");
  }

  ctx_free(col_nkmers);
  ctx_free(col_sum_covgs);

  // Close file (which zeros it)
  graph_file_close(&gfile);
  graph_header_dealloc(&hdr);

  return num_errors ? EXIT_FAILURE : EXIT_SUCCESS;
}
