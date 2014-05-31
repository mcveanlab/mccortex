#include "global.h"
#include "correct_aln_stats.h"
#include "util.h"

#define INIT_BUFLEN 1024

static inline void merge_arrays(uint64_t *restrict a,
                                const uint64_t *restrict b,
                                size_t n)
{
  size_t i;
  for(i = 0; i < n; i++) a[i] += b[i];
}

void correct_aln_stats_alloc(CorrectAlnStats *stats)
{
  // Gap Histogram
  stats->histgrm_len = INIT_BUFLEN;
  stats->gap_ins_histgrm = ctx_calloc(stats->histgrm_len, sizeof(uint64_t));
  stats->gap_err_histgrm = ctx_calloc(stats->histgrm_len, sizeof(uint64_t));
}

void correct_aln_stats_dealloc(CorrectAlnStats *stats)
{
  ctx_free(stats->gap_ins_histgrm);
  ctx_free(stats->gap_err_histgrm);
}

void correct_aln_stats_zero(CorrectAlnStats *stats)
{
  memset(stats->gap_err_histgrm, 0, stats->histgrm_len * sizeof(uint64_t));
  memset(stats->gap_ins_histgrm, 0, stats->histgrm_len * sizeof(uint64_t));
  stats->num_gap_attempts = 0;
  stats->num_gap_successes = 0;
  stats->num_paths_disagreed = 0;
  stats->num_gaps_too_short = 0;
}

void correct_aln_stats_merge(CorrectAlnStats *restrict dst,
                             CorrectAlnStats *restrict src)
{
  correct_aln_stats_cap(dst, src->histgrm_len);
  merge_arrays(dst->gap_err_histgrm, src->gap_err_histgrm, src->histgrm_len);
  merge_arrays(dst->gap_ins_histgrm, src->gap_ins_histgrm, src->histgrm_len);
  dst->num_gap_attempts += src->num_gap_attempts;
  dst->num_gap_successes += src->num_gap_successes;
  dst->num_paths_disagreed += src->num_paths_disagreed;
  dst->num_gaps_too_short += src->num_gaps_too_short;
}

void correct_aln_stats_cap(CorrectAlnStats *stats, size_t max_gap)
{
  if(stats->histgrm_len < max_gap) {
    max_gap = roundup2pow(max_gap);
    stats->gap_ins_histgrm = ctx_realloc(stats->gap_ins_histgrm, max_gap*sizeof(uint64_t));
    stats->gap_err_histgrm = ctx_realloc(stats->gap_err_histgrm, max_gap*sizeof(uint64_t));
    // Zero new memory
    size_t newmem = (max_gap-stats->histgrm_len)*sizeof(uint64_t);
    memset(stats->gap_ins_histgrm+stats->histgrm_len, 0, newmem);
    memset(stats->gap_err_histgrm+stats->histgrm_len, 0, newmem);
    stats->histgrm_len = max_gap;
  }
}

// Save gap size distribution
// insert_sizes is true if gaps are insert gaps,
//                 false if gaps are due to sequencing errors
void correct_aln_stats_dump(const char *path,
                            const uint64_t *arr, size_t arrlen,
                            size_t kmer_size, bool insert_sizes,
                            size_t nreads)
{
  ctx_assert(arrlen > 0);

  if(nreads == 0) { warn("No results to save."); return; }

  // Print summary statistics: min, mean, median, mode, max
  size_t i, min, max, total, ngaps = 0, mode = 0;
  max = total = arr[0];

  for(min = 0; min < arrlen && arr[min] == 0; min++) {}

  if(min == arrlen) {
    if(insert_sizes) status("No insert gaps traversed");
    else status("No seq error gaps traversed");
    return;
  }

  for(i = 1; i < arrlen; i++) {
    if(arr[i] > 0) max = i;
    if(arr[i] > arr[mode]) mode = i;
    ngaps += arr[i];
    total += arr[i] * i;
  }

  double mean = (double)total / ngaps;
  float median = find_hist_median(arr, arrlen, ngaps);

  size_t ninputs = insert_sizes ? nreads/2 : nreads;
  char ngaps_str[100], ninputs_str[100];
  ulong_to_str(ngaps, ngaps_str);
  ulong_to_str(ninputs, ninputs_str);

  status("%s size distribution: "
         "min: %zu mean: %.1f median: %.1f mode: %zu max: %zu",
         insert_sizes ? "Insert" : "Seq error gap",
         min, mean, median, mode, max);

  status("  Gaps per read%s: %s / %s [%.2f%%]",
         insert_sizes ? " pair" : "", ngaps_str, ninputs_str,
         (100.0*ngaps) / ninputs);

  FILE *fout;

  if((fout = fopen(path, "w")) == NULL) {
    warn("Cannot dump gapsize [cannot open: %s]", path);
    return;
  }

  fprintf(fout, "gap_in_kmers\tbp\tcount\n");

  if(arrlen > 0)
  {
    size_t start = 0, end = arrlen-1;

    while(start < arrlen && arr[start] == 0) start++;
    while(end > start && arr[end] == 0) end--;

    for(i = start; i <= end; i++) {
      fprintf(fout, "%4zu\t%4li\t%4zu\n",
              i, (long)i-(long)kmer_size, (size_t)arr[i]);
    }
  }

  status("Contig %s sizes dumped to %s\n",
         insert_sizes ? "insert" : "gap", path);

  fclose(fout);
}
