#include "global.h"
#include "loading_stats.h"
#include "util.h"

//
// Create/free/sum LoadingStats
//

void loading_stats_init(LoadingStats *stats)
{
  memset(stats, 0, sizeof(LoadingStats));
}

void loading_stats_merge(LoadingStats* dst, const LoadingStats* src)
{
  dst->num_files_loaded += src->num_files_loaded;
  dst->ctx_files_loaded += src->ctx_files_loaded;

  dst->num_se_reads += src->num_se_reads;
  dst->num_pe_reads += src->num_pe_reads;

  dst->num_good_reads += src->num_good_reads;
  dst->num_bad_reads += src->num_bad_reads;
  dst->num_dup_se_reads += src->num_dup_se_reads;
  dst->num_dup_pe_pairs += src->num_dup_pe_pairs;

  dst->total_bases_read += src->total_bases_read;
  dst->total_bases_loaded += src->total_bases_loaded;

  dst->contigs_loaded += src->contigs_loaded;
  dst->num_kmers_parsed += src->num_kmers_parsed;
  dst->num_kmers_loaded += src->num_kmers_loaded;
  dst->num_kmers_novel += src->num_kmers_novel;

  dst->num_of_colours_loaded += src->num_of_colours_loaded;
}

// @ht_num_kmers is the number of kmers loaded into the graph
void loading_stats_print_summary(const LoadingStats *stats, size_t ht_num_kmers)
{
  char se_num_str[100], pe_num_str[100], sepe_num_str[100];
  ulong_to_str(stats->num_se_reads, se_num_str);
  ulong_to_str(stats->num_pe_reads / 2, pe_num_str);
  ulong_to_str(stats->num_se_reads + stats->num_pe_reads, sepe_num_str);
  status("[SeqStats] single reads: %s; read pairs: %s; total: %s",
         se_num_str, pe_num_str, sepe_num_str);

  char nkmers_parsed_str[50], nkmers_loaded_str[50], nkmers_ht_str[50];
  ulong_to_str(stats->num_kmers_parsed, nkmers_parsed_str);
  ulong_to_str(stats->num_kmers_loaded, nkmers_loaded_str);
  ulong_to_str(ht_num_kmers,            nkmers_ht_str);

  // Estimate coverage
  if(ht_num_kmers > 0)
  {
    double covg      = (double)stats->num_kmers_parsed / ht_num_kmers;
    double covg_inf  = (double)stats->num_kmers_loaded / ht_num_kmers;
    status("[SeqStats] Input coverage: %.2fX (%s / %s) reconstructed: %.2fX (%s / %s)",
           covg,     nkmers_parsed_str, nkmers_ht_str,
           covg_inf, nkmers_loaded_str, nkmers_ht_str);
  }
  else {
    status("[SeqStats] Input kmers: %s reconstructed: %s",
           nkmers_parsed_str, nkmers_loaded_str);
  }

  if(stats->contigs_loaded > 0)
  {
    char klen_str[50];
    double mean_klen = (double)stats->num_kmers_loaded / stats->contigs_loaded;
    ulong_to_str((size_t)(mean_klen+0.5), klen_str);
    status("[SeqStats]  mean reconstructed contig length: %s (kmers)", klen_str);
  }
}
