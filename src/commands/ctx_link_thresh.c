#include "global.h"
#include "commands.h"
#include "cmd.h"
#include "util.h"
#include "file_util.h"
#include "clean_graph.h"

#define MAX_L_COUNT 2000

const char linkthresh_usage[] =
"usage: "CMD" linkthresh [options] <fdr> <dist> [in.csv]\n"
"\n"
"  Load a CSV of links and pick a threshold a given distance.\n"
"  CSV should have format <dist>,<count>\n"
"  <fdr> is the FDR threshold to use e.g. 0.001\n"
"\n"
"  -h, --help     This help message\n"
"  -q, --quiet    Silence status output normally printed to STDERR\n"
"  -z, --zero     Return zero if we cannot pick a threshold [default: throw error]\n"
"\n";

static struct option longopts[] =
{
  {"help", no_argument, NULL, 'h'},
  {"zero", no_argument, NULL, 'z'},
  {NULL, 0, NULL, 0}
};

// Returns two uint64. Expects format:
//     uint[,<tab><space>]uint\n
// will tolerate comment lines, empty lines, tab/sapce separated and a
// single header
static inline int parse_csv_line(const char *line, uint64_t *a, uint64_t *b)
{
  const char *sep = strchr(line, ',');
  char *endptr1, *endptr2;
  long aa, bb;
  if(sep == NULL) sep = strchr(line, '\t');
  if(sep == NULL) sep = strchr(line, ' ');
  if(sep == NULL) return -1;
  aa = strtoul(line, &endptr1, 10);
  bb = strtoul(sep+1, &endptr2, 10);
  if(endptr1 != sep || *endptr2 != '\0' || aa == 0 || bb == 0) return -1;
  *a = aa; *b = bb;
  return 0;
}

int ctx_linkthresh(int argc, char **argv)
{
  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  bool zero_on_error = false;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'z': cmd_check(!zero_on_error,cmd); zero_on_error = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("`"CMD" linkthresh -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  int rem_argc = argc - optind;
  char **rem_argv = argv + optind;

  if(rem_argc > 3 || rem_argc < 2) cmd_print_usage("Missing arguments");

  double fdr_limit = cmd_udouble_nonzero("<fdr>", rem_argv[0]);
  uint32_t dist = cmd_uint32_nonzero("<dist>", rem_argv[1]);
  if(fdr_limit >= 1) die("<fdr> should be 0<fdr<1 e.g. 0.001");
  const char *path = (rem_argc == 3 ? rem_argv[2] : "-");
  FILE *fin = futil_fopen(path, "r");

  uint64_t *counts = ctx_calloc(MAX_L_COUNT+1, sizeof(uint64_t));

  StrBuf line;
  strbuf_alloc(&line, 1024);
  bool exp_hdr = true;
  uint64_t a = 0, b = 0;

  while(strbuf_reset_readline(&line, fin))
  {
    strbuf_chomp(&line);
    if(line.b[0] != '#' && line.end > 0) {
      if(exp_hdr && (line.b[0] < '0' || line.b[0] > '9')) {}
      else {
        if(parse_csv_line(line.b, &a, &b) < 0)
          die("Bad CSV line: %s [%s]", line.b, path);
        if(a == dist) counts[MIN2(b, MAX_L_COUNT)]++;
      }
      exp_hdr = false;
    }
  }

  strbuf_dealloc(&line);

  double alpha = 0, beta = 0;
  int thresh = cleaning_pick_kmer_threshold(counts, MAX_L_COUNT+1, fdr_limit,
                                            &alpha, &beta);

  ctx_free(counts);
  fclose(fin);

  // Deal with error case where we cannot set threshold
  if(thresh < 0) {
    if(zero_on_error) thresh = 0;
    else {
      status("Cannot pick threshold at fdr = %f, dist = %u", fdr_limit, dist);
      return EXIT_FAILURE;
    }
  } else {
    status("alpha: %f beta: %f threshold: %i", alpha, beta, thresh);
  }

  printf("%i\n", thresh);
  return EXIT_SUCCESS;
}
