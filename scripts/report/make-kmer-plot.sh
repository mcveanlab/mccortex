#!/bin/bash
set -eou pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 <plot-covg-hist.R> <in.csv> <out.pdf>" 1>&2
  exit -1
fi

#in:  data/sample.kK.se.links.csv
#out: plots/sample.kK.se.links.pdf
script=$1
in=$2
out=$3

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
$script $in $out $CUTOFF $KCOV
