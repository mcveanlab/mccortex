#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "stream_buffer.h"

static void print_seqs(const char *seq, size_t len, const uint32_t arr[len][len])
{
  size_t i, j;
  for(i = 0; i < len; i++) {
    for(j = 0; j < i; j++) { // don't want to print matches twice
      if(arr[i][j] && (i+1==len || !arr[i+1][j+1])) {
        printf("%.*s\n", arr[i][j], seq+j+1-arr[i][j]);
      }
    }
  }
}

static void load_file(const char *path, char **buf, size_t *blen, size_t *bsize)
{
  FILE *fh = stdin;
  if(strcmp(path,"-") != 0 && (fh = fopen(path, "r")) == NULL) {
    fprintf(stderr, "Cannot read file: %s\n", path);
    exit(-1);
  }

  while(freadline(fh, buf, blen, bsize) > 0) {}

  if(strcmp(path,"-") != 0) fclose(fh);
}

int main(int argc, char **argv)
{
  if(argc != 2) {
    fprintf(stderr, "usage: %s <file>\n", argv[0]);
    return EXIT_FAILURE;
  }

  char *s = NULL;
  size_t len = 0, size = 0;
  load_file(argv[1], &s, &len, &size);

  size_t i, j, lensq = len * len;
  uint32_t *arr = calloc(lensq, sizeof(uint32_t));
  size_t idx, prev;

  for(i = 0; i < len; i++) arr[i*len] = (s[0] == s[i]);

  // i is row, j is column
  for(i = 1; i < len; i++) {
    idx = i*len+1;
    prev = (i-1)*len;
    for(j = 1; j < i; j++, idx++, prev++) {
      arr[idx] = s[i] == s[j] ? arr[prev]+1 : 0;
    }
  }

  // Print matches
  print_seqs(s, len, (const uint32_t (*)[len])arr);

  free(arr);
  free(s);

  return EXIT_SUCCESS;
}
