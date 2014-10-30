#ifndef CONTIG_CONFIDENCE_H_
#define CONTIG_CONFIDENCE_H_

#include "madcrowlib/madcrow_buffer.h"
madcrow_buffer(double_buf,DoubleBuffer,double);

typedef struct {
  DoubleBuffer table;
} ContigConfidenceTable;

// Call conf_table_destroy to release memory after calling this function
// void conf_table_load_csv(ContigConfidenceTable *conf_table,
//                          FILE *fh, const char *path);

// Call conf_table_destroy to release memory after calling this function
void conf_table_update_hist(ContigConfidenceTable *conf_table, size_t genome_size,
                            size_t *contig_hist, size_t hist_len);

// Call conf_table_destroy to release memory after calling this function
void conf_table_calc(ContigConfidenceTable *conf_table,
                     size_t max_read_len, double avg_bp_covg);

void conf_table_destroy(ContigConfidenceTable *conf_table);

double conf_table_lookup(const ContigConfidenceTable *table, size_t rlen);

void conf_table_print(const ContigConfidenceTable *table, FILE *fh);

void conf_table_save(const ContigConfidenceTable *table, const char *path);

#endif /* CONTIG_CONFIDENCE_H_ */
