#!/bin/bash

if [[ $# -lt 4 ]]; then
  echo "usage: $0 <outdir> <links|plain> <kmer ...>"
  exit
fi

out=$1
shift
strategy=$1
shift

# echo "out: $out/"
mkdir -p $out

for k in "$@"
do
  (echo k$k;
   dnacat -P k$k/perf.contigs.$strategy.rmdup.fa | \
     awk '{print length($0)}' | \
     sort -nr -;) > $out/$k.txt
done

a=`echo "{$@}" | tr ' ' ','`
x=`eval echo "$out/$a.txt"`
# echo $x
paste $x
