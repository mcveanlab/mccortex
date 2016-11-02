#!/bin/bash

set -eou pipefail

if [ $# -ne 3 ]; then
  >&2 echo "usage: $0 <fdr> <k> <tree.csv>" &&
  >&2 echo "  Pick threshold for cleaning links" && false
fi

fdr_limit=$1
k=$2
tree_csv=$3

maxk=$[ ( ($k + 31) / 32 ) * 32 - 1 ]
DIR=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )
CTX="$DIR/../../bin/mccortex $k"

thresh=$($CTX linkthresh -q --zero $fdr_limit $[$k+2] $tree_csv;
         $CTX linkthresh -q --zero $fdr_limit $[$k+3] $tree_csv;
         $CTX linkthresh -q --zero $fdr_limit $[$k+4] $tree_csv;
         $CTX linkthresh -q --zero $fdr_limit $[$k+5] $tree_csv;
         $CTX linkthresh -q --zero $fdr_limit $[$k+6] $tree_csv;)

# Print all 5 values
echo $thresh;
# Print median
echo $thresh | tr " " "\n" | sort -n | head -3 | tail -1
