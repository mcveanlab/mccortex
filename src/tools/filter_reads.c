

#include "global.h"

#include <zlib.h>

#include "string_buffer.h"
#include "seq_file.h"

#include "util.h"
#include "file_util.h"
#include "hash_table.h"
#include "binary_kmer.h"
#include "seq_reader.h"
#include "binary_format.h"

static const char usage[] =
"usage: filter_reads <in.ctx> <capacity> [FASTA|FASTQ] [OPTIONS]\n"
"  Filters reads based on which have a kmer in the graph. Options:\n"
"    --se_list <in.list> <out.fq.gz>\n"
"    --pe_list <pe.list1> <pe.list2> <out.1.fq.gz> <out.2.fq.gz>\n"
"  Can specify --se_list/--pe_list multiple times.\n"
"  Prints as gzipped FASTA unless FASTQ is specified\n";

typedef struct {
  char *in1, *in2;
  gzFile out1, out2;
  size_t num_of_reads_printed;
  void (*print)(const read_t *r, gzFile gz, int linewrap);
} AlignReadsData;

static hkey_t find_node(BinaryKmer bkmer, const dBGraph *db_graph)
{
  BinaryKmer tmpkey;
  db_node_get_key(bkmer, db_graph->kmer_size, tmpkey);
  return hash_table_find(&db_graph->ht, tmpkey);
}

static boolean read_touches_graph(const read_t *r, const dBGraph *db_graph,
                                  SeqLoadingStats *stats)
{
  boolean found = false;
  uint32_t kmer_size = db_graph->kmer_size, num_contigs = 0, kmers_loaded = 0;

  if(r->seq.end >= kmer_size)
  {
    size_t search_pos = 0, start, end = 0, i;
    BinaryKmer bkmer; Nucleotide nuc;

    while((start = seq_contig_start(r, search_pos, kmer_size, 0,0)) < r->seq.end &&
          !found)
    {
      end = seq_contig_end(r, start, kmer_size, 0, 0, &search_pos);
      stats->total_bases_loaded += end - start;
      num_contigs++;

      binary_kmer_from_str(r->seq.b + start, kmer_size, bkmer);
      kmers_loaded++;
      if(find_node(bkmer, db_graph) != HASH_NOT_FOUND) { found = true; break; }

      for(i = start+kmer_size; i < end; i++)
      {
        nuc = char_to_binary_nucleotide(r->seq.b[i]);
        binary_kmer_left_shift_one_base(bkmer, kmer_size);
        binary_kmer_set_last_nuc(bkmer, nuc);
        kmers_loaded++;
        if(find_node(bkmer, db_graph) != HASH_NOT_FOUND) { found = true; break; }
      }
    }
  }

  // Update stats
  stats->total_bases_read += r->seq.end;
  stats->kmers_loaded += kmers_loaded;
  stats->unique_kmers += kmers_loaded - found;
  if(num_contigs > 0) stats->total_good_reads++;
  else stats->total_bad_reads++;

  return found;
}

void filter_reads(read_t *r1, read_t *r2,
                  int qoffset1, int qoffset2,
                  SeqLoadingPrefs *prefs, SeqLoadingStats *stats, void *ptr)
{
  (void)qoffset1; (void)qoffset2;
  (void)prefs; (void)stats;

  AlignReadsData *data = (AlignReadsData*)ptr;
  const dBGraph *db_graph = prefs->db_graph;

  if(read_touches_graph(r1, db_graph, stats) ||
     (r2 != NULL && read_touches_graph(r2, db_graph, stats)))
  {
    if(r2 != NULL) {
      // Print paired-end
      gzFile out2 = data->out2 != NULL ? data->out2 : data->out1;
      data->print(r1, data->out1, 0);
      data->print(r2, out2, 0);
    }
    else {
      // Print single-ended
      data->print(r1, data->out1, 0);
    }
    data->num_of_reads_printed++;
  }
}

int main(int argc, char* argv[])
{
  if(argc < 3) print_usage(usage, NULL);

  unsigned long kmer_capacity;
  char *input_ctx_path = argv[1];

  if(!parse_entire_ulong(argv[2], &kmer_capacity))
    print_usage(usage, "Invalid kmer capacity: %lu\n", kmer_capacity);

  int argi, nextarg = 4;
  void (*print_func)(const read_t *r, gzFile gz, int linewrap) = seq_gzprint_fasta;

  if(strcasecmp(argv[nextarg], "FASTQ") == 0) {
    print_func = seq_gzprint_fastq;
    nextarg++;
  } else if(strcasecmp(argv[nextarg], "FASTA") == 0) {
    print_func = seq_gzprint_fasta;
    nextarg++;
  }

  // Test other input args
  // Check filelists are readable
  // Check output is writable

  for(argi = nextarg; argi < argc;)
  {
    if(strcasecmp(argv[argi], "--se_list") == 0)
    {
      if(argi + 2 >= argc) print_usage(usage, "Missing arguments");
      char *in = argv[argi+1], *out = argv[argi+2];
      if(!test_file_readable(in)) print_usage(usage, "Cannot read: %s", in);
      if(!test_file_writable(out)) print_usage(usage, "Cannot write: %s", out);
      argi += 3;
    }
    else if(strcasecmp(argv[argi], "--pe_list") == 0)
    {
      if(argi + 4 >= argc) print_usage(usage, "Missing arguments");
      char *in1 = argv[argi+1], *in2 = argv[argi+2];
      char *out1 = argv[argi+3], *out2 = argv[argi+4];
      if(!test_file_readable(in1)) print_usage(usage, "Cannot read: %s", in1);
      if(!test_file_readable(in2)) print_usage(usage, "Cannot read: %s", in2);
      if(!test_file_writable(out1)) print_usage(usage, "Cannot write: %s", out1);
      if(!test_file_writable(out2)) print_usage(usage, "Cannot write: %s", out2);
      argi += 5;
    }
    else print_usage(usage, "Unexpected argument: %s", argv[argi]);
  }

  // Probe binary to get kmer_size
  boolean is_binary = false;
  uint32_t kmer_size, num_of_cols;
  uint64_t num_kmers;

  if(!binary_probe(input_ctx_path, &is_binary, &kmer_size, &num_of_cols, &num_kmers))
    print_usage(usage, "Cannot read binary file: %s", input_ctx_path);
  else if(!is_binary)
    print_usage(usage, "Input binary file isn't valid: %s", input_ctx_path);

  // Load binary
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, kmer_capacity);

  binary_load(input_ctx_path, &db_graph,
              0, -1, true, false, NULL);

  // Set up parsing sequence
  SeqLoadingStats *stats = seq_loading_stats_create(0);
  SeqLoadingPrefs prefs = {.into_colour = 0, .load_seq = true,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false,
                           .load_binaries = false,
                           .update_ginfo = false,
                           .db_graph = &db_graph};

  //
  // Filtering reads
  //
  size_t total_reads_printed = 0;
  for(argi = nextarg; argi < argc;)
  {
    char is_se = (strcasecmp(argv[argi], "--se_list") == 0);
    char is_pe = (strcasecmp(argv[argi], "--pe_list") == 0);
    char *p1 = NULL, *p2 = NULL, *out1 = NULL, *out2 = NULL;
    size_t init_reads, reads_loaded;

    if(is_se) {
      p1 = argv[argi+1];
      out1 = argv[argi+2];
    }
    else if(is_pe) {
      p1 = argv[argi+1];
      p2 = argv[argi+2];
      out1 = argv[argi+3];
      out2 = argv[argi+4];
    }
    else die("Unexpected argument: %s", argv[argi]);

    AlignReadsData data = {p1, p2, NULL, NULL, 0, print_func};

    if((data.out1 = gzopen(out1, "w")) == NULL)
      die("Cannot write to: %s", out1);

    if(is_pe && (data.out2 = gzopen(out2, "w")) == NULL)
      die("Cannot write to: %s", out2);

    if(is_se) message("writing: %s\n", out1);
    else message("writing: %s, %s\n", out1, out2);

    init_reads = stats->total_good_reads + stats->total_bad_reads +
                 stats->total_dup_reads;

    parse_filelists(p1, p2, READ_FALIST,
                    &prefs, stats, filter_reads, &data);
    
    gzclose(data.out1);
    if(is_pe) gzclose(data.out2);

    total_reads_printed += data.num_of_reads_printed;
    reads_loaded = stats->total_good_reads + stats->total_bad_reads +
                   stats->total_dup_reads - init_reads;

    printf("  Printed %zu / %zu inputs\n",
           data.num_of_reads_printed, reads_loaded);

    if(is_se) argi += 3;
    else argi += 5;
  }

  size_t total_reads = stats->total_good_reads + stats->total_bad_reads +
                       stats->total_dup_reads;

  printf("Total printed %zu / %zu reads\n", total_reads_printed, total_reads);
  printf("Done\n");

  // free
  seq_loading_stats_free(stats);
  db_graph_dealloc(&db_graph);

  exit(EXIT_SUCCESS);
}
