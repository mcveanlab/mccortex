#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail

if [[ $# -ne 1 && $# -ne 2 ]];
then
  echo "usage ./cortex_to_ray.sh <input.ctx[:colour]>"
  exit -1;
fi

DIR=$( cd $( dirname ${BASH_SOURCE[0]} ) && cd .. && pwd )
CTX="$DIR/bin/mccortex31 view"

col=0
if [[ $# -eq 2 ]]
then
  col=$2

  # Check number of colours in binary
  bincols=`$CTX --info $1 | grep 'colours:' | grep -o '[0-9]*$'` || exit
  if [[ $col -ge $bincols ]]
  then
    echo "Binary only has $bincols colours (you requested $col)"
    exit -1
  fi
fi

$CTX --print_kmers $1 | awk 'BEGIN { col='$col' } {
  covg=$(2+col);
  edges=$(2+col+(NF-1)/2);
  x=substr(edges,0,4); y=substr(edges,5,8);
  gsub("\\.","",x); gsub("\\.","",y);
  print $1";"covg";"toupper(x)";"y
}'
