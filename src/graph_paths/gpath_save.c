#include "global.h"
#include "gpath_save.h"
#include "gpath_set.h"
#include "gpath_subset.h"
#include "binary_seq.h"
#include "cmd.h"
#include "util.h"

// #include <unistd.h> // gethostname()
#include <sys/utsname.h> // utsname()
#include <pwd.h>

// {
//   <JSON header>
// }
// <KMER> <num> .. (ignored)
// [FR] [nkmers] [njuncs] [nseen,nseen,nseen] [seq:ACAGT] .. (ignored)

// @subset is a temp variable that is reused each time
static void _gpath_save_node(hkey_t hkey, gzFile gzout,
                             GPathSubset *subset, const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;
  const size_t ncols = gpstore->gpset.ncols;
  GPath *first_gpath = gpath_store_fetch(gpstore, hkey);
  const GPath *gpath;
  size_t i;

  // Load and sort paths for given kmer
  gpath_subset_reset(subset);
  gpath_subset_load_llist(subset, first_gpath);
  gpath_subset_sort(subset);

  if(subset->list.len == 0) return;

  // Print "<kmer> <npaths>"
  BinaryKmer bkmer = db_graph->ht.table[hkey];
  char bkstr[MAX_KMER_SIZE+1];
  binary_kmer_to_str(bkmer, db_graph->kmer_size, bkstr);
  gzprintf(gzout, "%s %zu\n", bkstr, subset->list.len);

  char orchar[2] = {0};
  orchar[FORWARD] = 'F';
  orchar[REVERSE] = 'R';
  const uint8_t *nseenptr;
  size_t klen;

  for(i = 0; i < subset->list.len; i++)
  {
    gpath = subset->list.data[i];
    nseenptr = gpath_set_get_nseen(gpset, gpath);
    klen = gpath_set_get_klen(gpset, gpath);
    gzprintf(gzout, "%c %zu %u %u", orchar[gpath->orient], klen,
                                    gpath->num_juncs, (uint32_t)nseenptr[0]);

    for(i = 1; i < ncols; i++)
      gzprintf(gzout, ",%u", (uint32_t)nseenptr[i]);

    gzputc(gzout, ' ');
    binary_seq_gzprint(gpath->seq, gpath->num_juncs, gzout);
    gzputc(gzout, '\n');
  }
}

static void _gpath_save_hdr(gzFile gzout, const char *path,
                            cJSON **hdrs, size_t nhdrs,
                            const dBGraph *db_graph)
{
  const GPathStore *gpstore = &db_graph->gpstore;
  const GPathSet *gpset = &gpstore->gpset;

  // Construct cJSON
  cJSON *json = cJSON_CreateObject();

  char datestr[50];
  time_t date = time(NULL);
  strftime(datestr, sizeof(datestr), "%Y-%m-%d %H:%M:%S", localtime(&date));

  cJSON_AddStringToObject(json, "fileFormat", "ctp");
  cJSON_AddNumberToObject(json, "formatVersion", 2);
  cJSON_AddNumberToObject(json, "ncols", gpset->ncols);
  cJSON_AddNumberToObject(json, "kmer_size", db_graph->kmer_size);
  cJSON_AddNumberToObject(json, "num_kmers", gpstore->num_kmers_with_paths);
  cJSON_AddNumberToObject(json, "num_paths", gpstore->num_paths);
  cJSON_AddNumberToObject(json, "path_bytes", gpstore->path_bytes);
  cJSON *colours = cJSON_CreateArray();
  cJSON_AddItemToObject(json, "colours", colours);

  size_t i;
  for(i = 0; i < gpset->ncols; i++)
  {
    cJSON *sample = cJSON_CreateObject();
    cJSON_AddNumberToObject(sample, "colour", i);
    cJSON_AddStringToObject(sample, "sample", db_graph->ginfo[i].sample_name.buff);
    cJSON_AddItemToArray(colours, sample);
  }

  cJSON *commands = cJSON_CreateArray();
  cJSON_AddItemToObject(json, "commands", commands);

  // Add latest command
  char keystr[9];
  cJSON *command = cJSON_CreateObject();
  cJSON_AddItemToArray(commands, command);
  cJSON_AddStringToObject(command, "key", hex_rand_str(keystr, sizeof(keystr)));
  cJSON_AddStringToObject(command, "cmd", cmd_get_cmdline());
  cJSON_AddStringToObject(command, "cwd", cmd_get_cwd());

  // Get absolute path to output file
  char abspath[PATH_MAX + 1];
  if(realpath(path, abspath) != NULL)
    cJSON_AddStringToObject(command, "outpath", abspath);
  else
    cJSON_AddStringToObject(command, "outpath", path);

  cJSON_AddStringToObject(command, "date", datestr);
  cJSON_AddStringToObject(command, "cortex", CTX_VERSION);
  cJSON_AddStringToObject(command, "htslib", HTS_VERSION);
  cJSON_AddStringToObject(command, "zlib",   ZLIB_VERSION);

  // Get username
  // uid_t uid = geteuid();
  // struct passwd *pw = getpwuid(uid);
  // if(pw != NULL)
  //   cJSON_AddStringToObject(command, "user",     pw->pw_name);

  char username[1024+1];
  if(getlogin_r(username, sizeof(username)) == 0)
    cJSON_AddStringToObject(command, "user",     username);
  else
    warn("Couldn't get username");

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

  char *jstr = cJSON_Print(json);
  gzputs(gzout, jstr);
  gzputc(gzout, '\n');
  free(jstr);
  cJSON_Delete(json);
}

// @hdrs is array of JSON headers of input files
void gpath_save(gzFile gzout, const char *path,
                cJSON **hdrs, size_t nhdrs,
                dBGraph *db_graph)
{
  ctx_assert(gpath_set_has_nseen(&db_graph->gpstore.gpset));

  // Write header
  _gpath_save_hdr(gzout, path, hdrs, nhdrs, db_graph);

  // Print comments about the format
  gzputs(gzout, "\n");
  gzputs(gzout, "# This file was generated with Cortex\n");
  gzputs(gzout, "#   written by Isaac Turner\n");
  gzputs(gzout, "# \n");
  gzputs(gzout, "# Comment lines begin with a # and are ignored, but must come after the header\n");
  gzputs(gzout, "# Format is:\n");
  gzputs(gzout, "#   [kmer] [num_paths] ...(ignored)\n");
  gzputs(gzout, "#   [FR] [num_kmers] [num_juncs] [counts0,counts1,...] [juncs:ACAGT] ...(ignored)\n");
  gzputs(gzout, "\n");

  // Iterate over kmers writing paths
  GPathSubset subset;
  gpath_subset_alloc(&subset);
  gpath_subset_init(&subset, &db_graph->gpstore.gpset);
  HASH_ITERATE(&db_graph->ht, _gpath_save_node, gzout, &subset, db_graph);
  gpath_subset_dealloc(&subset);

  status("[GPathSave] Graph paths saved to %s", path);
}
