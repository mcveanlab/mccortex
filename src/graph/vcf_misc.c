#include "global.h"
#include "vcf_misc.h"
#include "util.h"
#include "file_util.h"
#include "json_hdr.h"

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

void vcf_hdrtxt_append_commands(cJSON *command, StrBuf *hdr, const char *path)
{
  bool first;
  for(; command != NULL; command = command->next)
  {
    cJSON *key  = json_hdr_get(command, "key",    cJSON_String,  path);
    cJSON *cmd  = json_hdr_get(command, "cmd",    cJSON_Array,   path);
    cJSON *cwd  = json_hdr_get(command, "cwd",    cJSON_String,  path);
    cJSON *prev = json_hdr_get(command, "prev",   cJSON_Array,   path);
    cJSON *ver  = json_hdr_try(command, "mccortex",cJSON_String, path);

    prev = prev->child; // result could be NULL
    if(prev && prev->type != cJSON_String) die("Invalid 'prev' field");
    strbuf_append_str(hdr, "##mccortex_");
    strbuf_append_str(hdr, key->valuestring);
    strbuf_append_str(hdr, "=<prev=\"");
    strbuf_append_str(hdr, prev ? prev->valuestring : "NULL");

    if(prev) {
      while((prev = prev->next) != NULL) {
        strbuf_append_str(hdr, ";");
        strbuf_append_str(hdr, prev->valuestring);
      }
    }
    strbuf_append_str(hdr, "\",cmd=\"");
    for(first = true, cmd = cmd->child; cmd; cmd = cmd->next, first = false) {
      if(!first) strbuf_append_char(hdr, ' ');
      strbuf_append_str(hdr, cmd->valuestring);
    }
    strbuf_append_str(hdr, "\",cwd=\"");
    strbuf_append_str(hdr, cwd->valuestring);
    strbuf_append_str(hdr, "\"");
    if(ver) {
      strbuf_append_str(hdr, ",version=\"");
      strbuf_append_str(hdr, ver->valuestring);
      strbuf_append_str(hdr, "\"");
    }
    strbuf_append_str(hdr, ">\n");
  }
}
