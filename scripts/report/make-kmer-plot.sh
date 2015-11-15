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

ROOT=`echo $in | awk '{gsub(/\.(raw|clean).cov.csv$/,"")}1'`

CUTOFF=`([[ -e $ROOT.kthresh ]] && cat $ROOT.kthresh) || echo 0`
KCOV=`([[ -e $ROOT.kmercov ]] && cat $ROOT.kmercov) || echo 0`

echo in=$in
echo out=$out
echo cutoff_file=$ROOT.kthresh
echo kcov_file=$ROOT.kmercov
echo CUTOFF=$CUTOFF
echo KCOV=$KCOV

set -o xtrace
$ctxdir/scripts/R/plot-covg-hist.R $in $out $CUTOFF $KCOV
