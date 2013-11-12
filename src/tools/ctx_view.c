#include "global.h"
#include <errno.h>

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "binary_kmer.h"
#include "hash_table.h"
#include "db_graph.h"
#include "db_node.h"
#include "graph_info.h"
#include "graph_format.h"
#include "graph_file_filter.h"

static const char usage[] =
"usage: "CMD" view [options] <in.ctx>\n"
" options:\n"
"   --kmers  Print kmers\n"
"   --check  Check kmers\n"
"   --info   Print info\n"
"\n"
" Default is [--info --check]\n";

static char* get_edges_str(Edges edges, char* kmer_colour_edge_str)
{
  int i;
  char str[] = "acgt";

  char left = edges >> 4;
  left = rev_nibble(left);
  char right = edges & 0xf;

  for(i = 0; i < 4; i++)
    kmer_colour_edge_str[i] = (left & (0x1 << i) ? str[i] : '.');

  for(i = 0; i < 4; i++)
    kmer_colour_edge_str[i+4] = toupper(right & (0x1 << i) ? str[i] : '.');

  kmer_colour_edge_str[8] = '\0';

  return kmer_colour_edge_str;
}

static void print_header(GraphFileHeader *h)
{
  printf("binary version: %u\n", h->version);
  printf("kmer size: %u\n", h->kmer_size);
  printf("bitfields: %u\n", h->num_of_bitfields);
  printf("colours: %u\n", h->num_of_cols);

  char num_kmers_str[50];
  ulong_to_str(h->num_of_kmers, num_kmers_str);
  printf("number of kmers: %s\n", num_kmers_str);
  printf("----\n");

  uint32_t i;
  for(i = 0; i < h->num_of_cols; i++)
  {
    GraphInfo *ginfo = h->ginfo + i;

    printf("Colour %i:\n", i);

    if(h->version >= 6)
    {
      // Version 6 only output
      printf("  sample name: '%s'\n", ginfo->sample_name.buff);
    }

    char total_sequence_str[100];
    ulong_to_str(ginfo->total_sequence, total_sequence_str);

    printf("  mean read length: %u\n", ginfo->mean_read_length);
    printf("  total sequence loaded: %s\n", total_sequence_str);

    if(h->version >= 6)
    {
      // Version 6 only output
      printf("  sequence error rate: %Lf\n", ginfo->seq_err);

      ErrorCleaning *ec = &ginfo->cleaning;
      printf("  tip clipping: %s\n", (ec->tip_clipping == 0 ? "no" : "yes"));

      printf("  remove low coverage supernodes: %s [threshold: %i]\n",
             ec->remv_low_cov_sups ? "yes" : "no",
             ec->remv_low_cov_sups_thresh);

      printf("  remove low coverage kmers: %s [threshold: %i]\n",
             ec->remv_low_cov_nodes ? "yes" : "no",
             ec->remv_low_cov_nodes_thresh);

      printf("  cleaned against graph: %s [against: '%s']\n",
             ec->is_graph_intersection ? "yes" : "no",
             ec->intersection_name.buff);
    }
  }
}

#define loading_warning(fmt,...) { num_warnings++; warn(fmt, ##__VA_ARGS__);}
#define loading_error(fmt,...) { num_errors++; warn(fmt, ##__VA_ARGS__);}

int ctx_view(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;
  if(argc == 0) print_usage(usage, NULL);

  boolean print_info = false, parse_kmers = false, print_kmers = false;
  size_t num_errors = 0, num_warnings = 0;
  int argi;

  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcmp(argv[argi],"--info") == 0) print_info = true;
    else if(strcmp(argv[argi],"--check") == 0) parse_kmers = true;
    else if(strcmp(argv[argi],"--kmers") == 0) print_kmers = true;
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  if(!print_info && !parse_kmers && !print_kmers)
    print_info = parse_kmers = true;

  if(argi+1 < argc) print_usage(usage, NULL);

  GraphFileReader file = INIT_GRAPH_READER;
  char *path = argv[argc-1];
  int ret = graph_file_open(&file, path, true);
  if(ret == 0) die("Cannot open file: %s", path);

  if(print_info)
  {
    char fsize_str[50];
    bytes_to_str(file.file_size, 0, fsize_str);
    printf("Loading file: %s\n", file.path.buff);
    printf("File size: %s\n", fsize_str);
    printf("----\n");
  }

  size_t ncols = graph_file_outncols(&file);

  GraphFileHeader outheader = INIT_GRAPH_FILE_HDR;
  graph_header_global_cpy(&outheader, &file.hdr);
  graph_header_alloc(&outheader, ncols);
  outheader.num_of_cols = ncols;

  size_t i, sum_covgs_read = 0, sum_seq_loaded = 0;
  size_t num_kmers_read = 0, num_all_zero_kmers = 0, num_zero_covg_kmers = 0;

  for(i = 0; i < file.ncols; i++) {
    graph_info_merge(outheader.ginfo + i, file.hdr.ginfo + file.cols[i]);
    sum_seq_loaded += outheader.ginfo[i].total_sequence;
  }

  // Print header
  if(print_info)
    print_header(&outheader);

  BinaryKmer bkmer;
  Covg covgs[ncols];
  Edges edges[ncols];
  char bkmerstr[MAX_KMER_SIZE+1], edgesstr[9];

  if(parse_kmers || print_kmers)
  {
    if(print_info && print_kmers) printf("----\n");

    for(; graph_file_read(&file, &bkmer, covgs, edges); num_kmers_read++)
    {
      // If kmer has no covg or edges -> don't load
      Covg keep_kmer = 0, covgs_sum = 0;
      for(i = 0; i < ncols; i++) {
        keep_kmer |= covgs[i] | edges[i];
        covgs_sum += covgs[i];
      }
      if(keep_kmer == 0) continue;

      sum_covgs_read += covgs_sum;

      /* Kmer Checks */
      // graph_file_read_kmer() already checks for:
      // 1. oversized kmers
      // 2. kmers with covg 0 in all colours

      // Check for all-zeros (i.e. all As kmer: AAAAAA)
      uint64_t kmer_words_or = 0;

      for(i = 0; i < file.hdr.num_of_bitfields; i++)
        kmer_words_or |= bkmer.b[i];

      if(kmer_words_or == 0)
      {
        if(num_all_zero_kmers == 1)
        {
          loading_error("more than one all 'A's kmers seen [index: %zu]\n",
                        num_kmers_read);
        }

        num_all_zero_kmers++;
      }

      // Check covg is 0 for all colours
      for(i = 0; i < ncols && covgs[i] == 0; i++);
      num_zero_covg_kmers += (i == ncols);

      // Print
      if(print_kmers)
      {
        binary_kmer_to_str(bkmer, file.hdr.kmer_size, bkmerstr);
        fputs(bkmerstr, stdout);

        // Print covgs
        for(i = 0; i < ncols; i++)
          fprintf(stdout, " %u", covgs[i]);

        // Print edges
        for(i = 0; i < ncols; i++) {
          fputc(' ', stdout);
          fputs(get_edges_str(edges[i], edgesstr), stdout);
        }

        fputc('\n', stdout);
      }
    }
  }

  // check for various reading errors
  if(errno != 0)
    loading_error("errno set [%i]: %s\n", (int)errno, strerror(errno));

  int err = ferror(file.fh);
  if(err != 0)
    loading_error("occurred after file reading [%i]\n", err);

  graph_file_close(&file);

  char num_str[50];

  if(print_kmers || parse_kmers)
  {
    if(num_kmers_read != file.hdr.num_of_kmers) {
      loading_warning("Expected %zu kmers, read %zu\n",
                      (size_t)file.hdr.num_of_kmers, (size_t)num_kmers_read);
    }

    if(num_all_zero_kmers > 1)
    {
      loading_error("%s all-zero-kmers seen\n",
                    ulong_to_str(num_all_zero_kmers, num_str));
    }

    if(num_zero_covg_kmers > 0)
    {
      loading_warning("%s kmers have no coverage in any colour\n",
                      ulong_to_str(num_zero_covg_kmers, num_str));
    }
  }

  if((print_kmers || parse_kmers) && print_info)
  {
    printf("----\n");
    printf("kmers read: %s\n", ulong_to_str(num_kmers_read, num_str));
    printf("covgs read: %s\n", ulong_to_str(sum_covgs_read, num_str));
    printf("seq loaded: %s\n", ulong_to_str(sum_seq_loaded, num_str));
  }

  if(print_info)
  {
    // Print memory stats
    uint64_t mem, capacity, num_buckets;
    uint8_t bucket_size;

    hash_table_cap(outheader.num_of_kmers / IDEAL_OCCUPANCY, true,
                   &num_buckets, &bucket_size);
    capacity = num_buckets * bucket_size;
    mem = capacity * (sizeof(BinaryKmer) +
          ncols * (sizeof(Covg) + sizeof(Edges)));

    char memstr[100], capacitystr[100], bucket_size_str[100], num_buckets_str[100];
    bytes_to_str(mem, 1, memstr);
    ulong_to_str(capacity, capacitystr);
    ulong_to_str(bucket_size, bucket_size_str);
    ulong_to_str(num_buckets, num_buckets_str);

    int mem_height = __builtin_ctzl((long)num_buckets);

    printf("----\n");
    printf("memory required: %s [capacity: %s]\n", memstr, capacitystr);
    printf("  bucket size: %s; number of buckets: %s\n",
            bucket_size_str, num_buckets_str);
    printf("  --kmer_size %u --mem_height %i --mem_width %i\n",
           file.hdr.kmer_size, mem_height, bucket_size);
  }

  if((print_kmers || parse_kmers) && print_info)
  {
    printf("----\n");
    if(num_warnings > 0 || num_errors > 0) {
      printf("Warnings: %zu; Errors: %zu\n",
              (size_t)num_warnings, (size_t)num_errors);
    }
    if(num_errors == 0)
      printf(num_warnings ? "Binary may be ok\n" : "Binary is valid\n");
  }

  graph_header_dealloc(&outheader);
  graph_file_dealloc(&file);

  return EXIT_SUCCESS;
}
