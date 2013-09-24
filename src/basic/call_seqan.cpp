#include <iostream>

#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include "call_seqan.h"

using namespace std;
using namespace seqan;

// was: match=0, mismatch=-1, gapextend=-1, gapopen=-2
// now: match=1, mismatch=-2, gapextend=-1, gapopen=-4
Score<int> scoringScheme(1, -2, -1, -4);
StringSet<DnaString> seq;
Graph<Alignment<StringSet<DnaString, Dependent<> > > > aliG(seq);
String<char> align;

void multiple_seq_align(char **seqs, int num, char *result, size_t *len)
{
  clear(seq);

  int i;
  for(i = 0; i < num; i++) {
    // cout << seqs[i] << endl;
    appendValue(seq, seqs[i]);
  }

  // Graph<Alignment<StringSet<DnaString, Dependent<> > > > aliG(seq);
  assignStringSet(aliG, seq);
  globalMsaAlignment(aliG, scoringScheme);
  convertAlignment(aliG, align);

  strcpy(result, toCString(align));
  *len = length(align);
}
