#!/bin/bash
set -eou pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <in.csv> <out.pdf>" 1>&2
  exit -1
fi

#in:  data/sample.kK.se.links.csv
#out: plots/sample.kK.se.links.pdf
ctxdir=$( cd $( dirname ${BASH_SOURCE[0]} ) && cd ../../ && pwd )
in=$1
out=$2

KMER=`echo "$in" | grep -oE 'k[0-9]+' | grep -oE '[0-9]+'`
CUTOFFFILE=`echo "$in" | awk '{gsub(/\.links\.csv$/,".links.thresh")}1'`
KCOVFILE=`echo "$in" | awk '{gsub(/\.(se|pe)\.links\.csv$/,".kmercov")}1'`
READLENFILE=`echo "$in" | awk '{gsub(/\.(se|pe)\.links\.csv$/,".readlen")}1'`

CUTOFF=`([[ -e $CUTOFFFILE ]] && cat $CUTOFFFILE) || echo 0`
KCOV=`([[ -e $KCOVFILE ]] && cat $KCOVFILE) || echo 0`
READLEN=`([[ -e $READLENFILE ]] && cat $READLENFILE) || echo 0`

echo KMER=$KMER
echo CUTOFFFILE=$CUTOFFFILE
echo KCOVFILE=$KCOVFILE
echo READLENFILE=$READLENFILE

set -o xtrace
$ctxdir/scripts/R/link-cov-heatmap.R $in $out $CUTOFF $KMER $KCOV $READLEN
