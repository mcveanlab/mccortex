#include "global.h"

#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "binary_format.h"

static const char usage[] =
"usage: ctx_build <kmer_size> <mem> <out.ctx> [OPTIONS]\n"
"  Build a cortex binary\n"
"\n"
"  Sequence options:\n"
"    --quality_score_threshold <qual>\n"
"      Filter for quality scores in the input file. [default: 0]\n"
"    --remove_pcr_duplicates\n"
"      Removes PCR duplicate reads by ignoring read pairs if both reads start at\n"
"      the same k-mer as a previous read\n"
"    --cut_homopolymers <bases>\n"
"      Breaks reads at homopolymers of length >= <bases>\n"
"      (i.e. max homopolymer in filtered read == threshold-1) [default: off]\n"
"\n"
"  Loading sequence:\n"
"    --se_list <se.list>               Load a list of FASTA/Q/BAM\n"
"    --pe_list <pe.list1> <pe.list2>   Load paired-end data\n"
"    --seq <in.fa|fq|sam> [in2.fq]     Load sequence data\n"
"      Consecutive sequence options are loaded into the same colour\n"
"\n"
"  Loading binaries\n"
"    --load_binary <data.colours>      Load a binary into new colour(s)\n"
"\n"
"  Loading multiple colours"
"    --colour_list <data.colours>      Load a colour list into new colour(s)\n";

int main(int argc, char **argv)
{
  
  return EXIT_SUCCESS;
}
