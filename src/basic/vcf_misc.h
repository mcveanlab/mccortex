#ifndef VCF_MISC_H_
#define VCF_MISC_H_

#include "htslib/vcf.h"

// VCF output type setting
// modes_htslib is for passing to hts_open, hsmodes_htslib is human readable
extern const char *modes_htslib[], *hsmodes_htslib[];

// reqtype if not null is used to set type
// reqtype: v=>vcf, z=>compressed vcf, b=>bcf, bu=>uncompressed bcf
// otherwise we guess from path file extension
// if both fail we return some default (uncompressed VCF)
// return index to address modes_htslib and hsmodes_htslib
int vcf_misc_get_outtype(const char *reqtype, const char *path);

void vcf_misc_hdr_add_cmd(bcf_hdr_t *hdr, const char *cmdline, const char *cwd);

// Find/add and then update a header record
void vcf_misc_add_update_hrec(bcf_hrec_t *hrec, char *key, char *val);

// Trim bases that match with the ref
static inline size_t trimmed_alt_lengths(const bcf1_t *v, size_t aid,
                                         size_t *rptr, size_t *aptr)
{
  const char *ref = v->d.allele[0], *alt = v->d.allele[aid];
  size_t rshift = 0;
  size_t reflen = strlen(ref);
  size_t altlen = strlen(alt);

  // Left trim
  while(reflen && altlen && *ref == *alt) {
    ref++; alt++; rshift++;
    reflen--; altlen--;
  }

  // Right trim
  while(reflen && altlen && ref[reflen-1] == alt[altlen-1]) {
    reflen--; altlen--;
  }

  *rptr = reflen;
  *aptr = altlen;

  return rshift;
}

#endif /* VCF_MISC_H_ */
