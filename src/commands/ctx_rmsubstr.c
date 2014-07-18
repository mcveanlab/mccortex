#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "seq_reader.h"
#include "kmer_occur.h"

const char rmsubstr_usage[] =
"usage: "CMD" rmsubstr [options] <in.fa> [in2.fq ...]\n"
"\n"
"  Remove duplicate sequences and those that occur as substrings of others.\n"
"\n"
"  -h, --help            This help message\n"
"  -q, --quiet           Silence status output normally printed to STDERR\n"
"  -f, --force           Overwrite output files\n"
"  -o, --out <out.txt>   Save output [default: STDOUT]\n"
"  -m, --memory <mem>    Memory to use\n"
"  -n, --nkmers <kmers>  Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>     Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -k, --kmer <kmer>     Kmer size must be odd ("QUOTE_VALUE(MAX_KMER_SIZE)" >= k >= "QUOTE_VALUE(MIN_KMER_SIZE)")\n"
"  -F, --fasta           Print output in FASTA format [default]\n"
"  -Q, --fastq           Print output in FASTQ format\n"
"  -P, --plain           Print output sequences one per line\n"
// "  -v, --invert                Print reads/read pairs with no kmer in graph\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"force",        no_argument,       NULL, 'f'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
// command specific
  {"kmer",         required_argument, NULL, 'k'},
  {"fasta",        no_argument,       NULL, 'F'},
  {"fastq",        no_argument,       NULL, 'Q'},
  {"plain",        no_argument,       NULL, 'P'},
  // {"invert",       no_argument,       NULL, 'v'},
  {NULL, 0, NULL, 0}
};

#if MAX_KMER_SIZE == 31
#  define DEFAULT_KMER 31
#else
#  define DEFAULT_KMER MIN_KMER_SIZE
#endif

#define OUTPUT_PLAIN 0
#define OUTPUT_FASTA 1
#define OUTPUT_FASTQ 2


// Returns true if a read is a substring of ANY read in the list or a complete
// match with a read before it in the list. Returns false otherwise.
//  1 => is substr
//  0 => not substr
// -1 => not enough bases of ACGT
static int _is_substr(const ReadBuffer *rbuf, size_t idx,
                      KOGraph kograph, const dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  const read_t *r = &rbuf->data[idx], *r2;
  const char *seq = r->seq.b, *seq2;
  size_t i, contig_start;

  contig_start = seq_contig_start(r, 0, kmer_size, 0, 0);
  if(contig_start >= r->seq.end) return -1;

  BinaryKmer bkmer = binary_kmer_from_str(seq+contig_start, kmer_size);
  dBNode node = db_graph_find(db_graph, bkmer);
  ctx_assert(node.key != HASH_NOT_FOUND);

  size_t num_hits = kograph_num(kograph, node.key);
  KOccur *hits = kograph_get(kograph, node.key);
  ctx_assert(num_hits > 0); // at least one hit (for this read!)

  for(i = 0; i < num_hits; i++)
  {
    if(hits[i].offset >= contig_start)
    {
      r2 = &rbuf->data[hits[i].chrom];
      seq2 = r2->seq.b + hits[i].offset - contig_start;

      // A read is a duplicate (i.e. return true) if it is a substring of ANY
      // read in the list or a complete match with a read before it in the list.
      // That is why we have: (hits[i].chrom < idx || r->seq.end != r2->seq.end)
      // since identical strings have equal length
      if(hits[i].chrom != idx && hits[i].orient == node.orient &&
         (hits[i].chrom < idx || r->seq.end != r2->seq.end) &&
         hits[i].offset + r->seq.end <= r2->seq.end &&
         strncasecmp(seq, seq2, r->seq.end) == 0)
      {
        return 1;
      }
    }
  }

  return 0;
}

int ctx_rmsubstr(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;
  size_t kmer_size = DEFAULT_KMER, num_of_threads = DEFAULT_NTHREADS;
  uint8_t output_format = OUTPUT_FASTA;
  size_t kmer_set = 0, output_format_set = 0;
  const char *output_file = NULL;

  // Arg parsing
  char cmd[100], shortopts[100];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'f': cmd_check(!futil_get_force(), cmd); futil_set_force(true); break;
      case 'o': cmd_check(!output_file, cmd); output_file = optarg; break;
      case 't': num_of_threads = cmd_uint32_nonzero(cmd, optarg); break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'k': kmer_set++; kmer_size = cmd_uint32(cmd, optarg); break;
      case 'F': output_format = OUTPUT_FASTA; output_format_set++; break;
      case 'Q': output_format = OUTPUT_FASTQ; output_format_set++; break;
      case 'P': output_format = OUTPUT_PLAIN; output_format_set++; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" rmsubstr -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(kmer_set > 1) cmd_print_usage("-k|--kmer specified more than once");
  if(output_format_set > 1) cmd_print_usage("Output format set more than once");

  if(optind >= argc)
    cmd_print_usage("Please specify at least one input graph file (.ctx)");

  size_t i, num_seq_files = argc - optind;
  char **seq_paths = argv + optind;
  seq_file_t **seq_files = ctx_calloc(num_seq_files, sizeof(seq_file_t*));
  size_t est_num_bases = 0;

  for(i = 0; i < num_seq_files; i++)
  {
    if((seq_files[i] = seq_open(seq_paths[i])) == NULL)
      die("Cannot read sequence file %s", seq_paths[i]);

    if(strcmp(seq_paths[i],"-") != 0) {
      off_t fsize = futil_get_file_size(seq_paths[i]);
      if(fsize < 0) warn("Cannot get file size: %s", seq_paths[i]);
      else {
        if(seq_is_fastq(seq_files[i]) || seq_is_sam(seq_files[i]))
          est_num_bases += fsize / 2;
        else
          est_num_bases += fsize;
      }
    }
  }

  // Use file sizes to decide on memory

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(BinaryKmer)*8 +
                  sizeof(KONodeList) + sizeof(KOccur) + // see kmer_occur.h
                  8; // 1 byte per kmer for each base to load sequence files

  kmers_in_hash = cmd_get_kmers_in_hash(memargs.mem_to_use,
                                        memargs.mem_to_use_set,
                                        memargs.num_kmers,
                                        memargs.num_kmers_set,
                                        bits_per_kmer,
                                        est_num_bases, est_num_bases,
                                        false, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  //
  // Open output file
  //
  if(output_file == NULL) output_file = "-";
  FILE *fout = futil_open_create(output_file, "w");

  //
  // Set up memory
  //
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, 1, 0, kmers_in_hash);
  db_graph.bktlocks = ctx_calloc(roundup_bits2bytes(db_graph.ht.num_of_buckets), 1);

  //
  // Load reference sequence into a read buffer
  //
  ReadBuffer rbuf;
  readbuf_alloc(&rbuf, 1024);
  seq_load_all_reads(seq_files, num_seq_files, &rbuf);

  // Check for reads too short
  for(i = 0; i < rbuf.len && rbuf.data[i].seq.end >= kmer_size; i++) {}
  if(i < rbuf.len)
    warn("Reads shorter than kmer size (%zu) will not be filtered", kmer_size);

  KOGraph kograph = kograph_create(rbuf.data, rbuf.len, true,
                                   num_of_threads, &db_graph);

  size_t num_reads = rbuf.len, num_reads_printed = 0, num_bad_reads = 0;

  // Loop over reads printing those that are not substrings
  int ret;
  for(i = 0; i < rbuf.len; i++) {
    ret = _is_substr(&rbuf, i, kograph, &db_graph);
    if(ret == 0) {
      switch(output_format) {
        case OUTPUT_PLAIN: fputs(rbuf.data[i].seq.b, fout); fputc('\n', fout); break;
        case OUTPUT_FASTA: seq_print_fasta(&rbuf.data[i], fout, 0); break;
        case OUTPUT_FASTQ: seq_print_fastq(&rbuf.data[i], fout, 0); break;
      }
      num_reads_printed++;
    }
    else if(ret == -1) num_bad_reads++;
  }

  char num_reads_str[100], num_reads_printed_str[100], num_bad_reads_str[100];
  ulong_to_str(num_reads, num_reads_str);
  ulong_to_str(num_reads_printed, num_reads_printed_str);
  ulong_to_str(num_bad_reads, num_bad_reads_str);

  status("Printed %s / %s (%.1f%%) to %s",
         num_reads_printed_str, num_reads_str,
         !num_reads ? 0.0 : (100.0 * num_reads_printed) / num_reads,
         futil_outpath_str(output_file));

  if(num_bad_reads > 0) {
    status("Bad reads: %s / %s (%.1f%%) - no kmer {ACGT} of length %zu",
           num_bad_reads_str, num_reads_str,
           (100.0 * num_bad_reads) / num_reads,
           kmer_size);
  }

  fclose(fout);
  kograph_free(kograph);

  // Free sequence memory
  for(i = 0; i < rbuf.len; i++) seq_read_dealloc(&rbuf.data[i]);
  readbuf_dealloc(&rbuf);
  ctx_free(seq_files);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
