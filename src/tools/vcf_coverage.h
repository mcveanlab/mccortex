#ifndef VCF_COVERAGE_H_
#define VCF_COVERAGE_H_

#include "db_graph.h"

#include "htslib/vcf.h"
#include "htslib/faidx.h"

#define DEFAULT_MAX_ALLELE_LEN 100
#define DEFAULT_MAX_GT_VARS 8

typedef struct {
  // Stats
  uint64_t nvcf_lines, nalts_read, nalts_loaded;
  uint64_t nalts_too_long, nalts_no_covg, nalts_with_covg;
  uint64_t ngt_kmers;
} VcfCovStats;

typedef struct {
  const char *kcov_ref_tag, *kcov_alt_tag;
  // Don't attempt to genotype alleles bigger than this
  // defaults to DEFAULT_MAX_ALLELE_LEN
  uint32_t max_allele_len;
  // 2^8 = 256 possible haplotypes
  // defaults to DEFAULT_MAX_GT_VARS
  uint32_t max_gt_vars;
  bool load_kmers_only;
} VcfCovPrefs;

void vcfcov_file(htsFile *vcffh, bcf_hdr_t *vcfhdr,
                 htsFile *outfh, bcf_hdr_t *outhdr,
                 const char *path, faidx_t *fai,
                 const size_t *samplehdrids,
                 const VcfCovPrefs *prefs,
                 VcfCovStats *stats,
                 dBGraph *db_graph);

#endif /* VCF_COVERAGE_H_ */
