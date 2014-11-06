#include "global.h"
#include "assemble_stats.h"

void assemble_contigs_stats_init(AssembleContigStats *stats)
{
  memset(stats, 0, sizeof(*stats));
  // Zero doubles
  stats->max_junc_density = 0.0;
  size_buf_capacity(&stats->lengths, 1<<20); // ~1 Million
  size_buf_capacity(&stats->junctns, 1<<20); // ~1 Million
}

void assemble_contigs_stats_destroy(AssembleContigStats *stats)
{
  size_buf_dealloc(&stats->lengths);
  size_buf_dealloc(&stats->junctns);
  memset(stats, 0, sizeof(*stats));
}

void assemble_contigs_stats_add(AssembleContigStats *stats,
                                const struct ContigStats *s)
{
  size_t i;

  // Update statistics
  for(i = 0; i < 2; i++) {
    stats->paths_held[MIN2(s->paths_held[i], AC_MAX_PATHS-1)]++;
    stats->paths_cntr[MIN2(s->paths_cntr[i], AC_MAX_PATHS-1)]++;
    stats->paths_held_max = MAX2(stats->paths_held_max, s->paths_held[i]);
    stats->paths_cntr_max = MAX2(stats->paths_cntr_max, s->paths_cntr[i]);
  }

  for(i = 0; i < GRPHWLK_NUM_STATES; i++)
    stats->grphwlk_steps[i] += s->wlk_steps[i];

  // Out degree
  stats->contigs_outdegree[s->outdegree_fw]++;
  stats->contigs_outdegree[s->outdegree_rv]++;

  size_buf_add(&stats->lengths, s->num_nodes);
  size_buf_add(&stats->junctns, s->njunc);

  stats->num_cycles += s->ncycles;
  stats->total_len  += s->num_nodes;
  stats->total_junc += s->njunc;

  stats->num_contigs_from_seed_paths += s->seed_path;
  stats->num_contigs_from_seed_kmers += s->seed_kmer;

  if(stats->num_contigs == 0) {
    stats->max_junc_density = (double)s->njunc / s->num_nodes;
  } else {
    stats->max_junc_density = MAX2(stats->max_junc_density,
                                   (double)s->njunc / s->num_nodes);
  }

  stats->num_contigs++;
}

void assemble_contigs_stats_merge(AssembleContigStats *dst,
                                  const AssembleContigStats *src)
{
  ctx_assert(dst->lengths.len == dst->junctns.len);
  ctx_assert(dst->lengths.len == dst->num_contigs);
  ctx_assert(src->lengths.len == src->junctns.len);
  ctx_assert(src->lengths.len == src->num_contigs);

  size_t i;

  size_buf_append(&dst->lengths, src->lengths.data, src->lengths.len);
  size_buf_append(&dst->junctns, src->junctns.data, src->junctns.len);

  dst->num_contigs += src->num_contigs;
  dst->num_cycles += src->num_cycles;
  dst->total_len += src->total_len;
  dst->total_junc += src->total_junc;

  for(i = 0; i < 5; i++)
    dst->contigs_outdegree[i] += src->contigs_outdegree[i];

  for(i = 0; i < AC_MAX_PATHS; i++) {
    dst->paths_held[i] += src->paths_held[i];
    dst->paths_cntr[i] += src->paths_cntr[i];
  }

  dst->paths_held_max = MAX2(dst->paths_held_max, src->paths_held_max);
  dst->paths_cntr_max = MAX2(dst->paths_cntr_max, src->paths_cntr_max);

  for(i = 0; i < GRPHWLK_NUM_STATES; i++)
    dst->grphwlk_steps[i] += src->grphwlk_steps[i];

  dst->max_junc_density = MAX2(dst->max_junc_density, src->max_junc_density);

  dst->num_contigs_from_seed_kmers += src->num_contigs_from_seed_kmers;
  dst->num_contigs_from_seed_paths += src->num_contigs_from_seed_paths;

  dst->num_reseed_abort += src->num_reseed_abort;
  dst->num_seeds_not_found += src->num_seeds_not_found;
}

#define PREFIX "[Assembled] "

static inline void pad_str(char *str, char c, size_t minlen)
{
  size_t len = strlen(str);
  if(len < minlen) {
    size_t shift = minlen - len;
    memmove(str+shift, str, len+1); // +1 to copy nul byte
    memset(str, c, shift);
  }
}

static inline void _print_grphwlk_state(const char *str, uint64_t nom,
                                        size_t denom)
{
  char nom_str[100], denom_str[100];
  ulong_to_str(nom, nom_str);
  ulong_to_str(denom, denom_str);
  pad_str(nom_str, ' ', 15);
  pad_str(denom_str, ' ', 15);
  status(PREFIX"  %s: %s / %s\t[ %2zu%% ]", str, nom_str, denom_str,
         !denom ? 0 : (size_t)(((100.0*nom)/denom)+0.5));
}

static inline void _print_path_dist(const uint64_t *hist, size_t n,
                                    const char *name, size_t num_contigs)
{
  char nout_str[100];
  size_t i;

  timestamp();
  message(PREFIX" %s: ", name);
  for(i = 0; i < n; i++) {
    message("\t%zu:%s [%zu%%]", i, ulong_to_str(hist[i], nout_str),
            (size_t)((100.0*hist[i])/(2.0*num_contigs)+0.5));
  }
  message("\n");
}


void assemble_contigs_stats_print(const AssembleContigStats *s)
{
  ctx_assert(s->lengths.len == s->junctns.len);
  ctx_assert(s->lengths.len == s->num_contigs);

  size_t i, ncontigs = s->num_contigs;

  if(ncontigs == 0) {
    status("[asm] No contigs assembled");
    return;
  }

  qsort(s->lengths.data, ncontigs, sizeof(s->lengths.data[0]), cmp_size);
  qsort(s->junctns.data, ncontigs, sizeof(s->junctns.data[0]), cmp_size);

  size_t len_n50, jnc_n50;
  size_t len_median, jnc_median, len_mean, jnc_mean;
  size_t len_min, len_max, jnc_min, jnc_max;

  // Calculate N50s
  len_n50 = calc_N50(s->lengths.data, ncontigs, s->total_len);
  jnc_n50 = calc_N50(s->junctns.data, ncontigs, s->total_junc);

  // Calculate medians, means
  len_median = MEDIAN(s->lengths.data, ncontigs);
  jnc_median = MEDIAN(s->junctns.data, ncontigs);
  len_mean = (double)s->total_len / ncontigs;
  jnc_mean = (double)s->total_junc / ncontigs;

  // Calculate min, max
  len_min = s->lengths.data[0];
  jnc_min = s->junctns.data[0];
  len_max = s->lengths.data[ncontigs-1];
  jnc_max = s->junctns.data[ncontigs-1];

  // Print number of contigs
  char num_contigs_str[50], reseed_str[50], seed_not_fnd_str[50];
  char seed_kmers_str[50], seed_paths_str[50];
  long_to_str(ncontigs, num_contigs_str);
  long_to_str(s->num_reseed_abort, reseed_str);
  long_to_str(s->num_seeds_not_found, seed_not_fnd_str);
  long_to_str(s->num_contigs_from_seed_kmers, seed_kmers_str);
  long_to_str(s->num_contigs_from_seed_paths, seed_paths_str);
  status(PREFIX"pulled out %s contigs, %s from seed kmers, %s from seed paths",
         num_contigs_str, seed_kmers_str, seed_paths_str);
  status(PREFIX"no-reseed aborted %s times", reseed_str);
  status(PREFIX"seed kmer not found %s times", seed_not_fnd_str);

  char len_min_str[50], len_max_str[50], len_total_str[50];
  char len_mean_str[50], len_median_str[50], len_n50_str[50];

  char jnc_min_str[50], jnc_max_str[50], jnc_total_str[50];
  char jnc_mean_str[50], jnc_median_str[50], jnc_n50_str[50];

  // Use ulong_to_str instead of num_to_str to get better accuracy
  // e.g. 966 instead of 1K
  ulong_to_str(len_mean, len_mean_str);
  ulong_to_str(jnc_mean, jnc_mean_str);
  ulong_to_str(len_median, len_median_str);
  ulong_to_str(jnc_median, jnc_median_str);
  ulong_to_str(len_n50, len_n50_str);
  ulong_to_str(jnc_n50, jnc_n50_str);
  ulong_to_str(len_min, len_min_str);
  ulong_to_str(jnc_min, jnc_min_str);
  ulong_to_str(len_max, len_max_str);
  ulong_to_str(jnc_max, jnc_max_str);
  ulong_to_str(s->total_len, len_total_str);
  ulong_to_str(s->total_junc, jnc_total_str);

  status(PREFIX"Lengths: mean: %s  median: %s  N50: %s  min: %s  max: %s  total: %s [kmers]",
         len_mean_str, len_median_str, len_n50_str, len_min_str, len_max_str, len_total_str);
  status(PREFIX"Junctions: mean: %s  median: %s  N50: %s  min: %s  max: %s  total: %s [out >1]",
         jnc_mean_str, jnc_median_str, jnc_n50_str, jnc_min_str, jnc_max_str, jnc_total_str);
  status(PREFIX"Max junction density: %.2f\n", s->max_junc_density);

  timestamp();
  message(PREFIX" Outdegree: ");
  char nout_str[50];

  for(i = 0; i <= 4; i++) {
    message("\t%zu:%s [%zu%%]", i, ulong_to_str(s->contigs_outdegree[i], nout_str),
            (size_t)((100.0*s->contigs_outdegree[i])/(2.0*ncontigs)+0.5));
  }
  message("\n");

  _print_path_dist(s->paths_held, AC_MAX_PATHS, "Paths held",    ncontigs);
  _print_path_dist(s->paths_cntr, AC_MAX_PATHS, "Paths counter", ncontigs);

  const uint64_t *states = s->grphwlk_steps;
  size_t nsteps = s->total_len - s->num_contigs, ncontigends = 2*s->num_contigs;
  status(PREFIX"Traversal succeeded because:");
  _print_grphwlk_state("Pop straight   ", states[GRPHWLK_POPFWD],       nsteps);
  _print_grphwlk_state("Col straight   ", states[GRPHWLK_COLFWD],       nsteps);
  _print_grphwlk_state("PopFork use col", states[GRPHWLK_POPFRK_COLFWD],nsteps);
  _print_grphwlk_state("Go paths       ", states[GRPHWLK_USEPATH],      nsteps);
  status(PREFIX"Traversal halted because:");
  _print_grphwlk_state("No coverage    ", states[GRPHWLK_NOCOVG],       ncontigends);
  _print_grphwlk_state("No colour covg ", states[GRPHWLK_NOCOLCOVG],    ncontigends);
  _print_grphwlk_state("No paths       ", states[GRPHWLK_NOPATHS],      ncontigends);
  _print_grphwlk_state("Paths split    ", states[GRPHWLK_SPLIT_PATHS],  ncontigends);
  _print_grphwlk_state("Missing paths  ", states[GRPHWLK_MISSING_PATHS],ncontigends);
  _print_grphwlk_state("Graph cycles   ", s->num_cycles,                ncontigends);

  size_t njunc = states[GRPHWLK_USEPATH] + states[GRPHWLK_NOPATHS] +
                 states[GRPHWLK_SPLIT_PATHS] + states[GRPHWLK_MISSING_PATHS];

  ctx_assert2(s->total_junc == states[GRPHWLK_USEPATH], "%zu vs %zu",
              (size_t)s->total_junc, (size_t)states[GRPHWLK_USEPATH]);

  status(PREFIX"Junctions:");
  _print_grphwlk_state("Paths resolved", states[GRPHWLK_USEPATH], njunc);
}
