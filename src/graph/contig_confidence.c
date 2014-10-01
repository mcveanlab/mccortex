#include "global.h"
#include "contig_confidence.h"

#include <math.h>

// expl(x) is the same as doing powl(M_E, x) a.k.a. 2**x for type long double
static double calc_confid(double bp_covg_depth, size_t read_length_bp,
                          size_t kmer_size)
{
  double lambda = bp_covg_depth / read_length_bp;
  double read_kmer_length = read_length_bp - kmer_size + 1;
  double power = (1 - expl(-lambda * read_kmer_length)) *
                 expl(-lambda * expl(-lambda * read_kmer_length));
  return power;
}

void conf_table_alloc(ContigConfidenceTable *conf_table,
                      size_t max_read_len, double avg_bp_covg)
{
  size_t i;

  status("Confidences for max. read length %zu and expected coverage %.2fX",
         max_read_len, avg_bp_covg);

  conf_table->read_len_max = max_read_len;
  conf_table->confids = ctx_calloc(max_read_len+1, sizeof(conf_table->confids[0]));
  conf_table->avg_bp_covg = avg_bp_covg;

  for(i = 0; i <= max_read_len; i++) {
    conf_table->confids[i] = calc_confid(avg_bp_covg, max_read_len, i);
  }
}

void conf_table_dealloc(ContigConfidenceTable *conf_table)
{
  ctx_free(conf_table->confids);
  memset(conf_table, 0, sizeof(ContigConfidenceTable));
}

double conf_table_lookup(const ContigConfidenceTable *table, size_t rlen)
{
  return rlen <= table->read_len_max ? table->confids[rlen] : 0;
}

void conf_table_print(const ContigConfidenceTable *table)
{
  size_t i;
  for(i = 0; i <= table->read_len_max; i++)
    printf("  %zu:%.2f\n", i, table->confids[i]);
  printf("\n");

  for(; i < table->read_len_max+10; i++)
    printf("  %zu:%.2f\n", i, calc_confid(table->avg_bp_covg, table->read_len_max, i));
}
