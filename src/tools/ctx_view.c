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
#include "binary_format.h"

static const char usage[] =
"usage: "CMD" view [options] <in.ctx>\n"
" options:\n"
"   --print_kmers\n"
"   --parse_kmers\n"
"   --print_info\n"
"\n"
" Default is [--print_info --parse_kmers]\n";

// A nibble is 4 bits (i.e. half a byte)
#define rev_nibble(x) (((x&0x1)<<3) | ((x&0x2)<<1) | ((x&0x4)>>1) | ((x&0x8)>>3))

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

static void print_header(BinaryFileHeader *h)
{
  message("binary version: %u\n", h->version);
  message("kmer size: %u\n", h->kmer_size);
  message("bitfields: %u\n", h->num_of_bitfields);
  message("colours: %u\n", h->num_of_cols);

  char num_kmers_str[50];
  ulong_to_str(h->num_of_kmers, num_kmers_str);
  message("number of kmers: %s\n", num_kmers_str);

  uint32_t i;
  for(i = 0; i < h->num_of_cols; i++)
  {
    GraphInfo *ginfo = h->ginfo + i;

    message("-- Colour %i --\n", i);

    if(h->version >= 6)
    {
      // Version 6 only output
      message("  sample name: '%s'\n", ginfo->sample_name.buff);
    }

    char total_sequence_str[100];
    ulong_to_str(ginfo->total_sequence, total_sequence_str);

    message("  mean read length: %u\n", ginfo->mean_read_length);
    message("  total sequence loaded: %s\n", total_sequence_str);
    
    if(h->version >= 6)
    {
      // Version 6 only output
      message("  sequence error rate: %Lf\n", ginfo->seq_err);

      ErrorCleaning *ec = &ginfo->cleaning;
      message("  tip clipping: %s\n", (ec->tip_clipping == 0 ? "no" : "yes"));

      message("  remove low coverage supernodes: %s [threshold: %i]\n",
              ec->remv_low_cov_sups ? "yes" : "no",
              ec->remv_low_cov_sups_thresh);

      message("  remove low coverage kmers: %s [threshold: %i]\n",
              ec->remv_low_cov_nodes ? "yes" : "no",
              ec->remv_low_cov_nodes_thresh);

      message("  cleaned against graph: %s [against: '%s']\n",
              ec->cleaned_against_another_graph ? "yes" : "no",
              ec->cleaned_against_graph_name.buff);
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
  char *in_ctx_path;
  int argi;

  for(argi = 0; argi < argc && argv[argi][0] == '-'; argi++)
  {
    if(strcmp(argv[argi],"--print_info") == 0) print_info = true;
    else if(strcmp(argv[argi],"--parse_kmers") == 0) parse_kmers = true;
    else if(strcmp(argv[argi],"--print_kmers") == 0) print_kmers = true;
    else print_usage(usage, "Unknown option: %s", argv[argi]);
  }

  if(!print_info && !parse_kmers && !print_kmers)
    print_info = parse_kmers = true;

  if(argi+1 < argc) print_usage(usage, NULL);

  in_ctx_path = argv[argc-1];

  FILE *in;

  char *split = strchr(in_ctx_path, ':');
  if(split != NULL) *split = '\0';

  // Get file size in case we need to count kmers
  off_t fsize = get_file_size(in_ctx_path);

  if((in = fopen(in_ctx_path, "r")) == NULL)
    die("Cannot open input path: %s", in_ctx_path);

  if(print_info)
  {
    char fsize_str[50];
    bytes_to_str(fsize, 0, fsize_str);
    message("Loading file: %s\n", in_ctx_path);
    message("File size: %s\n", fsize_str);
    message("----\n");
  }

  BinaryFileHeader outheader = {.capacity = 0}, inheader = {.capacity = 0};
  size_t hsize = binary_read_header(in, &inheader, in_ctx_path);

  uint32_t kmer_size = inheader.kmer_size;

  if(split != NULL) *split = ':';
  uint32_t num_of_cols = binary_get_num_colours(in_ctx_path, inheader.num_of_cols-1);
  uint32_t load_colours[num_of_cols];
  binary_parse_colour_array(in_ctx_path, load_colours, inheader.num_of_cols-1);

  binary_read_cpy_basic(&outheader, &inheader);
  outheader.num_of_cols = num_of_cols;
  binary_header_alloc(&outheader, num_of_cols);

  uint32_t i;
  uint64_t sum_of_covgs_read = 0, sum_of_seq_loaded = 0;

  for(i = 0; i < num_of_cols; i++) {
    graph_info_merge(outheader.ginfo + i, inheader.ginfo + load_colours[i]);
    sum_of_seq_loaded += outheader.ginfo[i].total_sequence;
  }

  if(outheader.version < 7)
  {
    size_t bytes_per_kmer = sizeof(BinaryKmer) +
                            outheader.num_of_cols * (sizeof(Covg) + sizeof(Edges));
    size_t bytes_remaining = fsize - hsize;
    outheader.num_of_kmers = (bytes_remaining / bytes_per_kmer);
    if(bytes_remaining % bytes_per_kmer != 0) {
      loading_warning("Truncated ctx binary: %s [bytes per kmer: %zu "
                      "remaining: %zu; fsize: %zu; header: %zu]",
                      in_ctx_path, bytes_per_kmer, bytes_remaining,
                      (size_t)fsize, hsize);
    }
  }

  // Print header
  if(print_info)
    print_header(&outheader);

  binary_header_dealloc(&inheader);
  binary_header_dealloc(&outheader);

  BinaryKmer bkmer;
  Covg kmercovgs[inheader.num_of_cols], covgs[num_of_cols];
  Edges kmeredges[inheader.num_of_cols], edges[num_of_cols];

  char bkmerstr[MAX_KMER_SIZE+1], edgesstr[9];

  uint64_t num_of_kmers_read = 0, num_of_oversized_kmers = 0,
           num_of_all_zero_kmers = 0, num_of_zero_covg_kmers = 0;

  // Check top word of each kmer
  int bits_in_top_word = 2 * (kmer_size % 32);
  uint64_t top_word_mask = (~(uint64_t)0) << bits_in_top_word;

  if(parse_kmers || print_kmers)
  {
    if(print_info && print_kmers) message("----\n");

    while(binary_read_kmer(in, &inheader, in_ctx_path, bkmer, covgs, edges))
    {
      // Collapse down colours
      if(split != NULL) {
        for(i = 0; i < num_of_cols; i++) {
          covgs[i] = kmercovgs[load_colours[i]];
          edges[i] = kmeredges[load_colours[i]];
        }
        // beware: if kmer has no covg or edges we still print it
      }

      for(i = 0; i < num_of_cols; i++)
        sum_of_covgs_read += covgs[i];

      /* Kmer Checks */
      // Check top bits of kmer
      if(bkmer[0] & top_word_mask)
      {
        if(num_of_oversized_kmers == 0) {
          loading_error("oversized kmer [index: %lu]\n",
                        (unsigned long)num_of_kmers_read);
        }

        num_of_oversized_kmers++;
      }

      // Check for all-zeros (i.e. all As kmer: AAAAAA)
      uint64_t kmer_words_or = 0;

      for(i = 0; i < outheader.num_of_bitfields; i++)
        kmer_words_or |= bkmer[i];

      if(kmer_words_or == 0)
      {
        if(num_of_all_zero_kmers == 1)
        {
          loading_error("more than one all 'A's kmers seen [index: %lu]\n",
                        (unsigned long)num_of_kmers_read);
        }

        num_of_all_zero_kmers++;
      }

      // Check covg is 0 for all colours
      for(i = 0; i < num_of_cols && covgs[i] == 0; i++);

      if(i == num_of_cols)
      {
        if(num_of_zero_covg_kmers == 0)
        {
          loading_warning("a kmer has zero coverage in all colours [index: %lu]\n",
                          (unsigned long)num_of_kmers_read);
        }

        num_of_zero_covg_kmers++;
      }

      // Print
      if(print_kmers)
      {
        binary_kmer_to_str(bkmer, kmer_size, bkmerstr);
        fputs(bkmerstr, stdout);

        // Print covgs
        for(i = 0; i < num_of_cols; i++)
          fprintf(stdout, " %u", covgs[i]);

        // Print edges
        for(i = 0; i < num_of_cols; i++) {
          fputc(' ', stdout);
          fputs(get_edges_str(edges[i], edgesstr), stdout);
        }

        fputc('\n', stdout);
      }

      num_of_kmers_read++;
    }
  }

  // check for various reading errors
  if(errno != 0)
    loading_error("errno set [%i]\n", (int)errno);

  int err = ferror(in);
  if(err != 0)
    loading_error("occurred after file reading [%i]\n", err);

  fclose(in);

  if(num_of_kmers_read != outheader.num_of_kmers) {
    loading_warning("Expected %zu kmers, read %zu\n",
                    (size_t)outheader.num_of_kmers, (size_t)num_of_kmers_read);
  }

  char num_str[50];

  if(num_of_all_zero_kmers > 1)
  {
    loading_error("%s all-zero-kmers seen\n",
                  ulong_to_str(num_of_all_zero_kmers, num_str));
  }

  if(num_of_oversized_kmers > 0)
  {
    loading_error("%s oversized kmers seen\n",
                  ulong_to_str(num_of_oversized_kmers, num_str));
  }

  if(num_of_zero_covg_kmers > 0)
  {
    loading_warning("%s kmers have no coverage in any colour\n",
                    ulong_to_str(num_of_zero_covg_kmers, num_str));
  }

  if((print_kmers || parse_kmers) && print_info)
  {
    message("----\n");
    message("kmers read: %s\n", ulong_to_str(num_of_kmers_read, num_str));
    message("covgs read: %s\n", ulong_to_str(sum_of_covgs_read, num_str));
    message("seq loaded: %s\n", ulong_to_str(sum_of_seq_loaded, num_str));
  }

  if((print_kmers || parse_kmers) && print_info)
  {
    message("----\n");
    if(num_warnings > 0 || num_errors > 0) {
      message("Warnings: %zu; Errors: %zu\n",
              (size_t)num_warnings, (size_t)num_errors);
    }
    if(num_errors == 0)
      message(num_warnings ? "Binary may be ok\n" : "Binary is valid\n");
  }

  return EXIT_SUCCESS;
}
