#!/bin/bash

set -eou pipefail
set -x

# usage: sga.sh <dir> <in.1.fq> <in.2.fq>


function abspath {
  echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CTXDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../.. && pwd )"

DNACAT=${CTXDIR}/libs/seq_file/bin/dnacat
PY_BREAK_VS_TRUTH=${CTXDIR}/scripts/python/break-contigs-vs-truth.py
REF=${CTXDIR}/results/data/chr22/chr22_28M_29M.fa

outdir=$1
IN1=`abspath "$2"`
IN2=`abspath "$3"`

mkdir -p $outdir
cd $outdir

for k in 21 31 41 51 61 71 81 91
do
  mkdir -p k$k
  cd k$k
  pwd
  $DIR/sga.sh $IN1 $IN2 $k
  $DNACAT -P $REF | $PY_BREAK_VS_TRUTH 21 assemble.m${k}-contigs.fa > stats.k${k}.out 2> stats.k${k}.txt
  cd ..
done
