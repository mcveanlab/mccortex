#include "global.h"
#include "loading_stats.h"

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
  dst->num_kmers_loaded += src->num_kmers_loaded;
  dst->num_kmers_novel += src->num_kmers_novel;

  dst->num_of_colours_loaded += src->num_of_colours_loaded;
}
