#ifndef CONTIG_CONFIDENCE_H_
#define CONTIG_CONFIDENCE_H_

typedef struct {
  size_t read_len_max;
  double *confids;
  double avg_bp_covg;
} ContigConfidenceTable;

void conf_table_alloc(ContigConfidenceTable *conf_table,
                      size_t max_read_len, double avg_bp_covg);

void conf_table_dealloc(ContigConfidenceTable *conf_table);

double conf_table_lookup(const ContigConfidenceTable *table, size_t rlen);

void conf_table_print(const ContigConfidenceTable *table);

#endif /* CONTIG_CONFIDENCE_H_ */
