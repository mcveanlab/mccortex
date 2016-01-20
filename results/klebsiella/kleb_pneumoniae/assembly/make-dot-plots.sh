#!/bin/bash

set -ou pipefail
set -o xtrace

#
# Make MUMMER dot plot
# Requires: mummer (nucmer, mummerplot), gnuplot, dnacat, ghostscript (ps2pdf)
#
# Clean up:
#   rm -rf links plain truthref

LINKCONTIGS=k31_links/klebPneu.rmdup.fa
PLAINCONTIGS=k31_plain/klebPneu.rmdup.fa
TRUTH=../truth/CAV1016.fa
REF=../ref/GCF_000016305.1_ASM1630v1_genomic.fa

CONTIGSTATS=~/mccortex/libs/bioinf-perl/fastn_scripts/contig_stats.pl
ASSEM_NCONTIGS=./get-max-covg.sh

REFLEN=$(dnacat -P $REF | awk '{x+=length($1)} END{print x}')

# Fix issues, re-run gnuplot
function fix_mummer_gp {
  local prefix=$1
  sed -i '.tmp' 's/^set mouse clipboardformat/\#set mouse clipboardformat/g' $prefix.gp
  gnuplot $prefix.gp
}

# make plot $prefix'_cov.pdf'
# make plot $prefix'_dot.pdf'
function mk_mummer_plots {
  local prefix=$1; shift
  # local otheropts=$2
  # coverage plot
  mummerplot -t postscript -c -p $prefix"_cov" $prefix.delta
  fix_mummer_gp $prefix"_cov"
  ps2pdf $prefix"_cov.ps" $prefix"_cov.pdf"
  # dot plot
  # removed: --breaklen 10000
  mummerplot -t postscript -l "$@" -p $prefix"_dot" $prefix.delta
  fix_mummer_gp $prefix"_dot"
  ps2pdf $prefix"_dot.ps" $prefix"_dot.pdf"
}

function get_longest_contigs {
  local seqin=$1
  local seqout=$2
  local reflen=$3
  local ncontigs=$($ASSEM_NCONTIGS $reflen $seqin | awk '{print $1}')
  set +o xtrace
  dnacat -P $seqin | awk '{print length($1),$1}' | sort -k1rn,2 | \
    awk '{print $2}' | head -$ncontigs | \
    dnacat -F -M <(i=0; while true; do echo "r$i"; ((i++)); done) > $seqout
  set -o xtrace
}

# mkdir -p links
# get_longest_contigs $LINKCONTIGS links/contigs.fa $REFLEN
# nucmer -mum -p links/links $TRUTH links/contigs.fa
# mk_mummer_plots links/links
# nucmer -mum -p links/linksswap links/contigs.fa $TRUTH
# mk_mummer_plots links/linksswap

# mkdir -p plain
# get_longest_contigs $PLAINCONTIGS plain/contigs.fa $REFLEN
# nucmer -mum -p plain/plain $TRUTH plain/contigs.fa
# mk_mummer_plots plain/plain
# nucmer -mum -p plain/plainswap plain/contigs.fa $TRUTH
# mk_mummer_plots plain/plainswap

# # Make plot between ref and truth
# mkdir -p truthref
# nucmer -mum -p truthref/truthref $TRUTH $REF
# mk_mummer_plots truthref/truthref
# nucmer -mum -p truthref/truthrefswap $REF $TRUTH
# mk_mummer_plots truthref/truthrefswap

$CONTIGSTATS links/contigs.fa > links/contig_stats.txt
$CONTIGSTATS plain/contigs.fa > plain/contig_stats.txt
