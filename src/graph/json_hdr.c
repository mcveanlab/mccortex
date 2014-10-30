#include "global.h"
#include "json_hdr.h"
#include "cmd.h"
#include "util.h"

// #include <unistd.h> // gethostname()
#include <sys/utsname.h> // utsname()
#include <pwd.h>


#define load_check(x,msg,...) if(!(x)) { die("[JSON] "msg, ##__VA_ARGS__); }

void json_hdr_read(FILE *fh, gzFile gz, const char *path, StrBuf *hdrstr)
{
  ctx_assert(fh == NULL || gz == NULL);

  size_t nread;
  strbuf_reset(hdrstr);

  nread = fh ? strbuf_readline(hdrstr, fh) : strbuf_gzreadline(hdrstr, gz);
  load_check(nread > 0, "Empty file: %s", path);
  load_check(hdrstr->b[0] == '{', "Expected JSON header: %s", path);

  size_t i, prev_offset = 0;
  size_t num_curly_open = 0, num_curly_close = 0; // '{' '}'
  size_t num_brkt_open = 0, num_brkt_close = 0; // '[' ']'
  bool in_string = false, escape_char = false; // '\'

  while(1)
  {
    for(i = prev_offset; i < hdrstr->end; i++) {
      if(in_string) {
        if(escape_char) escape_char = false;
        else if(hdrstr->b[i] == '\\') escape_char = true;
        else if(hdrstr->b[i] == '"') in_string = false;
      }
      else if(hdrstr->b[i] == '"') in_string = true;
      else if(hdrstr->b[i] == '{') num_curly_open++;
      else if(hdrstr->b[i] == '}') num_curly_close++;
      else if(hdrstr->b[i] == '[') num_brkt_open++;
      else if(hdrstr->b[i] == ']') num_brkt_close++;
    }
    prev_offset = hdrstr->end;

    if(num_curly_open == num_curly_close && num_brkt_open == num_brkt_close)
      break;

    // header is not finished yet - some checks
    load_check(num_curly_open > num_curly_close, "'}' before '{': %s", path);
    load_check(num_brkt_open >= num_brkt_close, "']' before '[': %s", path);
    load_check(hdrstr->end < MAX_JSON_HDR_BYTES, "Large JSON header: %s", path);

    // Read next line
    nread = fh ? strbuf_readline(hdrstr, fh) : strbuf_gzreadline(hdrstr, gz);
    load_check(nread > 0, "Premature end of JSON header: %s", path);
  }
}

// Add standard header fields to a json header
// Merge commands from input files @hdrs
// @param path is the path of the file we are writing to
void json_hdr_add_std(cJSON *json, const char *path,
                      cJSON **hdrs, size_t nhdrs,
                      const dBGraph *db_graph)
{
  // Add random id string
  #define FILE_ID_PREFIX "file:"
  #define FILE_ID_LEN 16
  char fileidstr[strlen(FILE_ID_PREFIX)+FILE_ID_LEN+1];
  strcpy(fileidstr, FILE_ID_PREFIX);
  hex_rand_str(fileidstr+strlen(FILE_ID_PREFIX), FILE_ID_LEN+1);
  cJSON_AddStringToObject(json, "file_id", fileidstr);

  cJSON *graph = cJSON_CreateObject();
  cJSON_AddItemToObject(json, "graph", graph);

  cJSON_AddNumberToObject(graph, "num_colours",        db_graph->num_of_cols);
  cJSON_AddNumberToObject(graph, "kmer_size",          db_graph->kmer_size);
  cJSON_AddNumberToObject(graph, "num_kmers_in_graph", db_graph->ht.num_kmers);

  cJSON *colours = cJSON_CreateArray();
  cJSON_AddItemToObject(graph, "colours", colours);

  size_t i;
  for(i = 0; i < db_graph->num_of_cols; i++)
  {
    bool cleaned_snodes = db_graph->ginfo[i].cleaning.cleaned_snodes;
    bool cleaned_tips   = db_graph->ginfo[i].cleaning.cleaned_tips;
    cJSON *sample = cJSON_CreateObject();
    cJSON_AddNumberToObject(sample, "colour", i);
    cJSON_AddStringToObject(sample, "sample", db_graph->ginfo[i].sample_name.b);
    cJSON_AddNumberToObject(sample, "total_sequence", db_graph->ginfo[i].total_sequence);
    cJSON_AddBoolToObject(sample, "cleaned_tips", cleaned_tips);
    if(cleaned_snodes) {
      cJSON_AddNumberToObject(sample, "cleaned_supernodes",
                              db_graph->ginfo[i].cleaning.clean_snodes_thresh);
    } else {
      cJSON_AddBoolToObject(sample, "cleaned_supernodes", false);
    }
    cJSON_AddItemToArray(colours, sample);
  }

  cJSON *commands = cJSON_CreateArray();
  cJSON_AddItemToObject(json, "commands", commands);

  // Add latest command
  char keystr[4+8+1];
  strcpy(keystr, "cmd:");
  hex_rand_str(keystr+4, 9);

  cJSON *command = cJSON_CreateObject();
  cJSON_AddItemToArray(commands, command);
  cJSON_AddStringToObject(command, "key", keystr);

  // Add command line arguments
  const char **argv = cmd_get_argv();
  int argc = cmd_get_argc();
  cJSON *cmdargs = cJSON_CreateStringArray(argv, argc);
  cJSON_AddItemToObject(command, "cmd", cmdargs);

  cJSON_AddStringToObject(command, "cwd", cmd_get_cwd());

  // Get absolute path to output file
  char abspath[PATH_MAX + 1];
  if(realpath(path, abspath) != NULL)
    cJSON_AddStringToObject(command, "outpath", abspath);
  else {
    status("Warning: Cannot get absolute path for: %s", path);
    cJSON_AddStringToObject(command, "outpath", path);
  }

  char datestr[50];
  time_t date = time(NULL);
  strftime(datestr, sizeof(datestr), "%Y-%m-%d %H:%M:%S", localtime(&date));

  cJSON_AddStringToObject(command, "date", datestr);
  cJSON_AddStringToObject(command, "cortex", CTX_VERSION);
  cJSON_AddStringToObject(command, "htslib", HTS_VERSION);
  cJSON_AddStringToObject(command, "zlib",   ZLIB_VERSION);

  // Get username
  // struct passwd *pw = getpwuid(geteuid());
  // if(pw != NULL)
  //   cJSON_AddStringToObject(command, "user",     pw->pw_name);
  struct passwd pwd;
  struct passwd *pwd_result;
  char username[1024+1] = "noname";

  // getpwuid may succeed when getlogin fails (e.g. when running not in terminal)
  if(getlogin_r(username, sizeof(username)) == 0 ||
     (getpwuid_r(geteuid(), &pwd, username, sizeof(username), &pwd_result) == 0 &&
      pwd_result != NULL))
  {
    cJSON_AddStringToObject(command, "user", username);
  } else {
    status("Warning: Cannot get username for JSON header");
  }

  // Get system info
  struct utsname sysdata;
  if(uname(&sysdata) != -1) {
    cJSON_AddStringToObject(command, "host",      sysdata.nodename);
    cJSON_AddStringToObject(command, "os",        sysdata.sysname);
    cJSON_AddStringToObject(command, "osrelease", sysdata.release);
    cJSON_AddStringToObject(command, "osversion", sysdata.version);
    cJSON_AddStringToObject(command, "hardware",  sysdata.machine);
  }

  // char hostname[2048];
  // if(gethostname(hostname, sizeof(hostname)) != -1)
  //   cJSON_AddStringToObject(command, "host", hostname);

  cJSON *prev_list = cJSON_CreateArray();
  cJSON_AddItemToObject(command, "prev", prev_list);

  size_t j;
  CharPtrBuffer cmdbuf;
  char_ptr_buf_alloc(&cmdbuf, nhdrs*8);

  // for each file loaded, link to first command
  for(i = 0; i < nhdrs; i++) {
    cJSON *first = cJSON_GetObjectItem(hdrs[i], "commands");
    if(first != NULL && first->type == cJSON_Array && first->child != NULL) {
      cJSON *prev_command = first->child;

      // Add link to first prev_command
      cJSON *key = cJSON_GetObjectItem(prev_command, "key");
      if(key != NULL && key->type == cJSON_String)
        cJSON_AddItemToArray(prev_list, cJSON_CreateString(key->valuestring));

      // Then add their commands to the list
      for(; prev_command != NULL; prev_command = prev_command->next)
      {
        // Only add if "key" not duplicate
        key = cJSON_GetObjectItem(prev_command, "key");
        if(key != NULL && key->type == cJSON_String)
        {
          for(j = 0; j < cmdbuf.len; j++)
            if(strcmp(cmdbuf.data[j], key->valuestring) == 0)
              break;

          if(j == cmdbuf.len) {
            cJSON_AddItemToArray(commands, cJSON_Duplicate(prev_command, 1));
            char_ptr_buf_add(&cmdbuf, key->valuestring);
          }
        }
      }
    }
  }

  char_ptr_buf_dealloc(&cmdbuf);
}

void json_hdr_gzprint(cJSON *json, gzFile gzout)
{
  char *jstr = cJSON_Print(json);
  gzputs(gzout, jstr);
  gzputc(gzout, '\n');
  free(jstr);
}

cJSON* json_hdr_get(cJSON *json, const char *field, int type, const char *path)
{
  cJSON *obj = cJSON_GetObjectItem(json, field);
  if(obj == NULL)
    die("Cannot find field in JSON header: %s [path: %s]", field, path);
  if(obj->type != type)
    die("JSON field not of correct type: %s [path: %s]", field, path);
  return obj;
}

long json_hdr_demand_int(cJSON *json, const char *field, const char *path)
{
  cJSON *obj = json_hdr_get(json, field, cJSON_Number, path);
  return obj->valueint;
}

size_t json_hdr_demand_uint(cJSON *json, const char *field, const char *path)
{
  cJSON *obj = json_hdr_get(json, field, cJSON_Number, path);
  if(obj->valueint < 0) {
    die("JSON field should not be less than zero: %s (%li) [path: %s]",
        field, obj->valueint, path);
  }
  return obj->valueint;
}

size_t json_hdr_get_kmer_size(cJSON *json, const char *path)
{
  cJSON *graph = json_hdr_get_graph(json, path);

  long val = json_hdr_demand_int(graph, "kmer_size", path);
  if(val < MIN_KMER_SIZE || val > MAX_KMER_SIZE || !(val & 1)) {
    die("kmer size is not an odd int between %i..%i: %li [%s]",
        MIN_KMER_SIZE, MAX_KMER_SIZE, val, path);
  }
  return val;
}

size_t json_hdr_get_ncols(cJSON *json, const char *path)
{
  cJSON *graph = json_hdr_get_graph(json, path);
  size_t val = json_hdr_demand_uint(graph, "num_colours", path);
  if(val < 1 || val > 100000) die("Invalid number of colours: %zu", val);
  return val;
}

cJSON* json_hdr_get_curr_cmd(cJSON *json, const char *path)
{
  cJSON *cmds = json_hdr_get(json, "commands", cJSON_Array, path);
  if(cmds->child == NULL) die("No 'commands' field in header");
  return cmds->child;
}
