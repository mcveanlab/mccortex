#include "global.h"
#include "vcf_misc.h"
#include "util.h"
#include "file_util.h"

// v=>vcf, z=>compressed vcf, b=>bcf, bu=>uncompressed bcf
const char *modes_htslib[] = {"wv","wz","wbu","wb"};
const char *hsmodes_htslib[] = {"uncompressed VCF", "compressed VCF",
                                "uncompressed BCF", "compressed BCF"};

// reqtype: v=>vcf, z=>compressed vcf, b=>bcf, bu=>uncompressed bcf
int vcf_misc_get_outtype(const char *reqtype, const char *path)
{
  int mode = -1;
  if(reqtype) {
    const char *ot = reqtype;
    if(!strcmp(ot,"vcf") || !strcmp(ot,"uvcf") || !strcmp(ot,"v")) mode = 0;
    else if(!strcmp(ot,"vcfgz") || !strcmp(ot,"z")) mode = 1;
    else if(!strcmp(ot,"ubcf")  || !strcmp(ot,"u")) mode = 2;
    else if(!strcmp(ot,"bcf")   || !strcmp(ot,"b")) mode = 3;
    else die("Unkown output type: %s", ot);
  }
  else if(path) {
    if(     futil_path_has_extension(path,".vcf"))    mode = 0;
    else if(futil_path_has_extension(path,".vcfgz"))  mode = 1;
    else if(futil_path_has_extension(path,".vcf.gz")) mode = 1;
    else if(futil_path_has_extension(path,".ubcf"))   mode = 2;
    else if(futil_path_has_extension(path,".bcf"))    mode = 3;
  }
  // default to uncompressed VCF
  return mode < 0 ? 0 : mode;
}

void vcf_misc_hdr_add_cmd(bcf_hdr_t *hdr, const char *cmdline, const char *cwd)
{
  char keystr[8], timestr[100];
  time_t tnow;
  time(&tnow);
  strftime(timestr, sizeof(timestr), "%Y%m%d-%H:%M:%S", localtime(&tnow));
  StrBuf sbuf;
  strbuf_alloc(&sbuf, 1024);
  strbuf_sprintf(&sbuf, "##mccortex_%s=<prev=\"NULL\",cmd=\"%s\",cwd=\"%s\","
                        "datetime=\"%s\",version="CTX_VERSION">\n",
                 hex_rand_str(keystr, sizeof(keystr)),
                 cmdline, cwd, timestr);
  bcf_hdr_append(hdr, sbuf.b);
  strbuf_dealloc(&sbuf);
}

// Find/add and then update a header record
void vcf_misc_add_update_hrec(bcf_hrec_t *hrec, char *key, char *val)
{
  int keyidx = bcf_hrec_find_key(hrec, key);
  if(keyidx < 0) {
    status("Adding sample key [%s] => [%s]", key, val);
    bcf_hrec_add_key(hrec, key, strlen(key));
    keyidx = hrec->nkeys-1;
  }
  bcf_hrec_set_val(hrec, keyidx, val, strlen(val), 0); // 0 => not quoted
}
