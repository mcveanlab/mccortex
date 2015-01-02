#ifndef VAR_ALIGNMENT_H_
#define VAR_ALIGNMENT_H_

// Returns msa_len if no more variants
size_t var_aln_msa_start(const char **alleles, size_t num_alleles,
                         size_t msa_len, size_t pos);

size_t var_aln_msa_end(const char **alleles, size_t num_alleles,
                       size_t msa_len, size_t pos);

// Copy alignment from aln to allele, removing '-'
void var_aln_strip_allele(StrBuf *allele, const char *aln, size_t len);

void var_aln_left_right_shift(const char *fl5p, size_t fl5plen,
                              const char *allele, size_t allelelen,
                              const char *fl3p, size_t fl3plen,
                              size_t *left, size_t *right);

#endif /* VAR_ALIGNMENT_H_ */
