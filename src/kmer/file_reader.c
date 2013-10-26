#include "global.h"
#include "file_reader.h"

//
// Create/free/sum SeqLoadingStats
//

SeqLoadingStats* seq_loading_stats_create(size_t readlen_arrsize)
{
  SeqLoadingStats *stats = calloc2(1, sizeof(SeqLoadingStats));
  stats->readlen_count_array_size = readlen_arrsize;

  if(readlen_arrsize > 0)
    stats->readlen_count_array = calloc2(readlen_arrsize, sizeof(unsigned long));
  else
    stats->readlen_count_array = NULL;

  return stats;
}

void seq_loading_stats_sum(SeqLoadingStats* dst, SeqLoadingStats* src)
{
  dst->se_colourlists_loaded += src->se_colourlists_loaded;
  dst->pe_colourlist_pairs_loaded += src->pe_colourlist_pairs_loaded;
  dst->se_filelists_loaded += src->se_filelists_loaded;
  dst->pe_filelist_pairs_loaded += src->pe_filelist_pairs_loaded;
  dst->num_files_loaded += src->num_files_loaded;
  dst->binaries_loaded += src->binaries_loaded;

  dst->num_se_reads += src->num_se_reads;
  dst->num_pe_reads += src->num_pe_reads;
  dst->total_good_reads += src->total_good_reads;
  dst->total_bad_reads += src->total_bad_reads;
  dst->total_dup_reads += src->total_dup_reads;
  dst->total_bases_read += src->total_bases_read;
  dst->total_bases_loaded += src->total_bases_loaded;
  dst->kmers_loaded += src->kmers_loaded;
  dst->unique_kmers += src->unique_kmers;
  dst->contigs_loaded += src->contigs_loaded;

  // Used for binaries and colourlists
  dst->num_of_colours_loaded += src->num_of_colours_loaded;

  uint64_t i;
  uint64_t limit = MIN2(dst->readlen_count_array_size,
                        src->readlen_count_array_size);

  for(i = 0; i < limit; i++)
  {
    dst->readlen_count_array[i] += src->readlen_count_array[i];
  }

  if(dst->readlen_count_array_size < src->readlen_count_array_size)
  {
    uint64_t sum = 0;

    for(i = dst->readlen_count_array_size; i < src->readlen_count_array_size; i++)
    {
      sum += src->readlen_count_array[i];
    }

    dst->readlen_count_array[dst->readlen_count_array_size-1] += sum;
  }
}

void seq_loading_stats_free(SeqLoadingStats* stats)
{
  if(stats->readlen_count_array != NULL) free(stats->readlen_count_array);
  free(stats);
}
