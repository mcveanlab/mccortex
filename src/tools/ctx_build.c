#include "global.h"

#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_info.h"
#include "graph_format.h"
#include "file_reader.h"
#include "seq_reader.h"

static const char usage[] =
"usage: "CMD" build [options] <out.ctx>\n"
"  Build a cortex graph.  \n"
"\n"
"  Options:\n"
"    -m <mem>               Memory to use (e.g. 100G or 12M)\n"
"    -k <kmer>              Kmer size\n"
"    --sample <name>        Sample name (required before any seq args)\n"
"    --seq <in.fa|fq|sam>   Load sequence data\n"
"    --seq2 <in1> <in2>     Load paired end sequence data\n"
"    --fq_threshold <qual>  Filter quality scores [default: 0 (off)]\n"
"    --fq_offset <offset>   FASTQ ASCII offset    [default: 0 (auto-detect)]\n"
"    --cut_hp <len>         Breaks reads at homopolymers >= <len> [default: off]\n"
"    --remove_pcr           Remove (or keep) PCR duplicate reads [default: keep]\n"
"    --keep_pcr             Don't do PCR duplicate removal\n"
"    --load_graph <in.ctx>  Load samples from a graph file (.ctx)\n"
"\n"
"  PCR duplicate removal works by ignoring read pairs (PE-only) if both reads\n"
"  start at the same k-mer as any previous read. Carried out per sample, not \n"
"  per file.\n"
"\n"
"  --sample <name> is required before sequence input can be loaded.\n"
"  Consecutive sequence options are loaded into the same colour.\n"
"\n"
"  --load_graph argument can have colours specifed e.g. in.ctx:0,6-8 will load\n"
"  samples 0,6,7,8.  Graphs are loaded into new colours.\n"
"\n"
" Example:\n"
"   "CMD" build -k 31 -m 60G --fq_threshold 5 --cut_hp 8 \\\n"
"               --sample bob --seq bob.bam \\\n"
"               --sample dylan --remove_pcr --seq dylan.sam \\\n"
"                              --keep_pcr --fq_offset 33 \\\n"
"                              --seq2 dylan.1.fq.gz dylan.2.fq.gz \\\n"
"               bob.dylan.k31.ctx\n"
"\n"
"  See `"CMD" join` to combine .ctx files\n";

static void update_ginfo(GraphInfo *ginfo, SeqLoadingStats *stats,
                         uint64_t *bases_loaded, uint64_t *contigs_loaded)
{
  uint64_t bases = stats->total_bases_loaded - *bases_loaded;
  uint64_t contigs = stats->contigs_loaded - *contigs_loaded;
  graph_info_update_contigs(ginfo, bases, contigs);
  *bases_loaded = stats->total_bases_loaded;
  *contigs_loaded = stats->contigs_loaded;
}

static void print_prefs(const SeqLoadingPrefs *prefs)
{
  char minQual[30] = "off", fqOffset[30] = "off", hpCutoff[30] = "auto-detect";

  if(prefs->quality_cutoff > 0)
    sprintf(minQual, "%u", prefs->quality_cutoff);
  if(prefs->ascii_fq_offset > 0)
    sprintf(fqOffset, "%u", prefs->ascii_fq_offset);
  if(prefs->homopolymer_cutoff > 0)
    sprintf(hpCutoff, "%u", prefs->homopolymer_cutoff);

  status("[prefs] FASTQ minQual: %s; FASTQ offset: %s; cut homopolymers: %s; "
         "remove PCR duplicates SE: %s, PE: %s\n", minQual, fqOffset, hpCutoff,
         prefs->remove_dups_se ? "yes" : "no",
         prefs->remove_dups_pe ? "yes" : "no");
}

int ctx_build(CmdArgs *args)
{
  cmd_accept_options(args, "mnk", usage);
  cmd_require_options(args, "k", usage);
  int argc = args->argc;
  char **argv = args->argv;
  if(argc < 3) print_usage(usage, NULL);

  size_t kmer_size = args->kmer_size;

  const char *out_path = argv[argc-1];
  size_t output_colours = 0;
  boolean sample_named = false, sample_used = false, remove_pcr_used = false;

  // Validate arguments
  int argi, argend = argc-1;
  uint32_t tmp;
  GraphFileReader ctxfile = INIT_GRAPH_READER;
  seq_file_t *seqfiles[argc];
  size_t num_sf = 0, sf = 0;

  for(argi = 0; argi < argend; argi++)
  {
    if(strcmp(argv[argi],"--fq_threshold") == 0) {
      if(argi + 1 >= argend)
        print_usage(usage, "--fq_threshold <qual> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &tmp) || tmp > 128)
        die("Invalid --fq_threshold argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_offset") == 0) {
      if(argi + 1 >= argend)
        print_usage(usage, "--fq_offset <offset> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &tmp) || tmp > 128)
        die("Invalid --fq_offset argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--cut_hp") == 0) {
      if(argi + 1 >= argend)
        print_usage(usage, "--cut_hp <len> requires an argument");
      if(!parse_entire_uint(argv[argi+1], &tmp))
        die("Invalid --cut_hp argument: %s", argv[argi+1]);
      argi += 1;
    }
    else if(!strcmp(argv[argi],"--remove_pcr")) { remove_pcr_used = true; }
    else if(!strcmp(argv[argi],"--keep_pcr")) {}
    else if(strcmp(argv[argi],"--seq") == 0) {
      if(!sample_named)
        print_usage(usage, "Please use --sample <name> before giving sequence");
      if(argi + 1 >= argend)
        print_usage(usage, "--seq <file> requires an argument");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read --seq file: %s", argv[argi+1]);
      argi += 1;
      sample_used = true;
    }
    else if(strcmp(argv[argi],"--seq2") == 0) {
      if(!sample_named)
        print_usage(usage, "Please use --sample <name> before giving sequence");
      if(argi + 2 >= argend)
        print_usage(usage, "--seq2 <file1> <file2> requires two arguments");
      if((seqfiles[num_sf++] = seq_open(argv[argi+1])) == NULL)
        die("Cannot read first --seq2 file: %s", argv[argi+1]);
      if((seqfiles[num_sf++] = seq_open(argv[argi+2])) == NULL)
        die("Cannot read second --seq2 file: %s", argv[argi+2]);
      argi += 2;
      sample_used = true;
    }
    else if(strcmp(argv[argi],"--load_graph") == 0) {
      if(argi + 1 >= argend) print_usage(usage, "--load_graph requires an arg");

      char *path = argv[argi+1];
      int ret = graph_file_open(&ctxfile, path, false);

      if(ret == 0)
        print_usage(usage, "Cannot read input graph file: %s", path);
      else if(ret < 0)
        print_usage(usage, "Input graph file isn't valid: %s", path);

      if(ctxfile.hdr.kmer_size != kmer_size) {
        print_usage(usage, "Input graph kmer_size doesn't match [%u vs %zu]",
                    ctxfile.hdr.kmer_size, kmer_size);
      }
      argi++;
      output_colours += graph_file_usedcols(&ctxfile);
      sample_named = false;
    }
    else if(strcmp(argv[argi],"--sample") == 0) {
      if(argi + 1 >= argend)
        print_usage(usage, "--sample <name> requires an argument");
      if(output_colours > 0 && !sample_used)
        warn("Empty colour (maybe you intended this)");
      if(!strcmp(argv[argi+1],"undefined") || !strcmp(argv[argi+1],"noname"))
        die("--sample %s is not a good name!", argv[argi+1]);
      argi++;
      output_colours++;
      sample_named = true;
      sample_used = false;
    }
    else {
      print_usage(usage, "Unknown command: %s", argv[argi]);
    }
  }

  if(!futil_is_file_writable(out_path))
    die("Cannot write to file: %s", out_path);

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash;

  bits_per_kmer = ((sizeof(Covg) + sizeof(Edges)) * 8 + remove_pcr_used*2) *
                  output_colours;
  kmers_in_hash = cmd_get_kmers_in_hash(args, bits_per_kmer, 0, true);

  status("Writing %zu colour graph to %s\n", output_colours, out_path);

  // Create db_graph
  dBGraph db_graph;
  db_graph_alloc(&db_graph, kmer_size, output_colours, output_colours, kmers_in_hash);
  db_graph.col_edges = calloc2(db_graph.ht.capacity * output_colours, sizeof(Edges));
  db_graph.col_covgs = calloc2(db_graph.ht.capacity * output_colours, sizeof(Covg));

  size_t kmer_words = round_bits_to_words64(db_graph.ht.capacity);

  if(remove_pcr_used)
    db_graph.readstrt = calloc2(kmer_words*2, sizeof(uint64_t));

  hash_table_print_stats(&db_graph.ht);

  // Parse arguments, load
  SeqLoadingStats *stats = seq_loading_stats_create(1000);
  SeqLoadingPrefs prefs = {.db_graph = &db_graph, .into_colour = -1,
                           .quality_cutoff = 0, .ascii_fq_offset = 0,
                           .homopolymer_cutoff = 0,
                           .remove_dups_se = false, .remove_dups_pe = false};

  read_t r1, r2;
  if(seq_read_alloc(&r1) == NULL || seq_read_alloc(&r2) == NULL)
    die("Out of memory");

  uint64_t bases_loaded = 0, contigs_loaded = 0;
  boolean show_prefs = true;

  for(argi = 0; argi < argend; argi++)
  {
    if(strcmp(argv[argi],"--sample") == 0) {
      prefs.into_colour++;
      status("[sample] %zu: %s\n", prefs.into_colour, argv[argi+1]);
      strbuf_set(&db_graph.ginfo[prefs.into_colour].sample_name, argv[argi+1]);
      if(db_graph.readstrt != NULL && prefs.into_colour > 0)
        memset(db_graph.readstrt, 0, 2 * kmer_words * sizeof(uint64_t));
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_threshold") == 0) {
      parse_entire_uint(argv[argi+1], &tmp);
      prefs.quality_cutoff = tmp;
      show_prefs = true;
      argi += 1;
    }
    else if(strcmp(argv[argi],"--fq_offset") == 0) {
      parse_entire_uint(argv[argi+1], &tmp);
      prefs.ascii_fq_offset = tmp;
      show_prefs = true;
      argi += 1;
    }
    else if(strcmp(argv[argi],"--remove_pcr") == 0) {
      prefs.remove_dups_pe = true;
      show_prefs = true;
    }
    else if(strcmp(argv[argi],"--keep_pcr") == 0) {
      prefs.remove_dups_pe = false;
      show_prefs = true;
    }
    else if(strcmp(argv[argi],"--cut_hp") == 0) {
      parse_entire_uint(argv[argi+1], &tmp);
      prefs.homopolymer_cutoff = tmp;
      show_prefs = true;
      argi += 1;
    }
    else if(strcmp(argv[argi],"--seq") == 0) {
      if(show_prefs) print_prefs(&prefs);
      show_prefs = false;
      seq_parse_se_sf(seqfiles[sf++], &r1, &r2, &prefs, stats,
                      seq_load_into_db_graph, NULL);
      update_ginfo(&db_graph.ginfo[prefs.into_colour],
                   stats, &bases_loaded, &contigs_loaded);
      argi += 1;
    }
    else if(strcmp(argv[argi],"--seq2") == 0) {
      if(show_prefs) print_prefs(&prefs);
      show_prefs = false;
      seq_parse_pe_sf(seqfiles[sf], seqfiles[sf+1], &r1, &r2, &prefs, stats,
                      seq_load_into_db_graph, NULL);
      update_ginfo(&db_graph.ginfo[prefs.into_colour],
                   stats, &bases_loaded, &contigs_loaded);
      argi += 2;
      sf += 2;
    }
    else if(strcmp(argv[argi],"--load_graph") == 0) {
      graph_file_open(&ctxfile, argv[argi+1], true);
      graph_load(&ctxfile, &prefs, stats);
      prefs.into_colour += graph_file_usedcols(&ctxfile);
      graph_file_close(&ctxfile);
      argi += 1;
    }
    else {
      die("Unknown command: %s", argv[argi]);
    }
  }

  seq_read_dealloc(&r1);
  seq_read_dealloc(&r2);
  seq_loading_stats_free(stats);

  hash_table_print_stats(&db_graph.ht);

  status("Dumping graph...\n");
  graph_file_save_mkhdr(out_path, &db_graph, CTX_GRAPH_FILEFORMAT, NULL,
                        0, output_colours);

  free(db_graph.col_covgs);
  free(db_graph.col_edges);
  if(db_graph.readstrt != NULL) free(db_graph.readstrt);

  db_graph_dealloc(&db_graph);
  graph_file_dealloc(&ctxfile);

  return EXIT_SUCCESS;
}
