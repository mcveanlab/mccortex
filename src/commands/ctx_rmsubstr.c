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
"  Options:\n"
"    -m <mem>    Memory to use\n"
"    -n <kmers>  Hash size\n"
"    -k <kmer>   Output file [default: STDOUT]\n"
"    --fasta     Print output in FASTA format [default]\n"
"    --fastq     Print output in FASTQ format\n"
"    --plain     Print output sequences one per line\n";

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
// Returns false if the read is shorter than kmer_size.
static bool _is_substr(const ReadBuffer *rbuf, size_t idx,
                       KOGraph kograph, const dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  const read_t *r = &rbuf->data[idx], *r2;
  const char *seq = r->seq.b, *seq2;

  if(r->seq.end < kmer_size) return false;

  size_t i;
  BinaryKmer bkmer = binary_kmer_from_str(seq, kmer_size);
  dBNode node = db_graph_find(db_graph, bkmer);
  ctx_assert(node.key != HASH_NOT_FOUND);

  size_t num_hits = kograph_num(kograph, node.key);
  KOccur *hits = kograph_get(kograph, node.key);
  ctx_assert(num_hits > 0); // at least one hit (for this read!)

  for(i = 0; i < num_hits; i++)
  {
    r2 = &rbuf->data[hits[i].chrom];
    seq2 = r2->seq.b + hits[i].offset;

    // A read is a duplicate (i.e. return true) if it is a substring of ANY
    // read in the list or a complete match with a read before it in the list.
    // That is why we have: (hits[i].chrom < idx || r->seq.end != r2->seq.end)
    // since identical strings have equal length
    if(hits[i].chrom != idx && hits[i].orient == node.orient &&
       (hits[i].chrom < idx || r->seq.end != r2->seq.end) &&
       hits[i].offset + r->seq.end <= r2->seq.end &&
       strncasecmp(seq, seq2, r->seq.end) == 0)
    {
      return true;
    }
  }

  return false;
}

int ctx_rmsubstr(CmdArgs *args)
{
  int argc = args->argc;
  char **argv = args->argv;

  size_t kmer_size = DEFAULT_KMER;
  uint8_t output_format = OUTPUT_FASTA;
  size_t kmer_set = 0, output_format_set = 0;

  int argi;
  for(argi = 0; argi < argc && argv[argi][0] == '-' && argv[argi][1]; argi++)
  {
    if(!strcasecmp(argv[argi], "--kmer") | !strcasecmp(argv[argi], "-k"))
    {
      if(argi+1 >= argc || !parse_entire_size(argv[argi+1], &kmer_size) || 
         (kmer_size & 1) == 0)
      {
        cmd_print_usage("%s <kmer> requires an odd integer value %i >= K >= %i",
                        argv[argi], MIN_KMER_SIZE, MAX_KMER_SIZE);
      }
      kmer_set = true;
      argi++;
    }
    else if(strcasecmp(argv[argi], "--plain") == 0) {
      output_format = OUTPUT_PLAIN;
      output_format_set++;
    }
    else if(strcasecmp(argv[argi], "--fasta") == 0) {
      output_format = OUTPUT_FASTA;
      output_format_set++;
    }
    else if(strcasecmp(argv[argi], "--fastq") == 0) {
      output_format = OUTPUT_FASTQ;
      output_format_set++;
    }
    else cmd_print_usage("Unknown option: %s", argv[argi]);
  }

  if(kmer_set > 1) cmd_print_usage("--kmer|-k specified more than once");
  if(output_format_set > 1) cmd_print_usage("Output format set more than once");

  if(argi >= argc)
    cmd_print_usage("Please specify at least one input graph file (.ctx)");

  size_t i, num_seq_files = argc - argi;
  char **seq_paths = argv + argi;
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
  size_t num_of_threads = args->max_work_threads;

  bits_per_kmer = sizeof(KONodeList) + sizeof(KOccur) + // see kmer_occur.h
                  8; // 1 byte per kmer for each base to load sequence files

  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer,
                                        est_num_bases, est_num_bases,
                                        false, &graph_mem);

  cmd_check_mem_limit(args, graph_mem);

  //
  // Open output file
  //
  const char *output_path = args->output_file_set ? args->output_file : "-";
  FILE *fout = futil_open_output(output_path);

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

  for(i = 0; i < num_seq_files && rbuf.data[i].seq.end >= kmer_size; i++) {}
  if(i < num_seq_files)
    warn("Reads shorter than kmer size (%zu) will not be filtered", kmer_size);

  KOGraph kograph = kograph_create(rbuf.data, rbuf.len, true,
                                   num_of_threads, &db_graph);

  size_t num_reads_start = rbuf.len, num_reads_end = 0;

  // Loop over reads printing those that are not substrings
  for(i = 0; i < rbuf.len; i++) {
    if(!_is_substr(&rbuf, i, kograph, &db_graph)) {
      switch(output_format) {
        case OUTPUT_PLAIN: fputs(rbuf.data[i].seq.b, fout); fputc('\n', fout); break;
        case OUTPUT_FASTA: seq_print_fasta(&rbuf.data[i], fout, 0); break;
        case OUTPUT_FASTQ: seq_print_fastq(&rbuf.data[i], fout, 0); break;
      }
      num_reads_end++;
    }
  }

  char num_reads_start_str[100], num_reads_end_str[100];
  ulong_to_str(num_reads_start, num_reads_start_str);
  ulong_to_str(num_reads_end, num_reads_end_str);

  status("Printed %s / %s (%.1f%%) to %s",
         num_reads_end_str, num_reads_start_str,
         (100.0 * num_reads_end) / num_reads_start,
         futil_outpath_str(output_path));

  fclose(fout);
  kograph_free(kograph);
  ctx_free(db_graph.bktlocks);
  db_graph_dealloc(&db_graph);

  // Free sequence memory
  for(i = 0; i < rbuf.len; i++) seq_read_dealloc(&rbuf.data[i]);
  readbuf_dealloc(&rbuf);
  ctx_free(seq_files);

  return EXIT_SUCCESS;
}
