#!/bin/bash

set -euo pipefail

function usage {
  echo "usage: $0 <chr22.fa>" >&2
  echo "  Compare buffered vs unbuffered read times" >&2
  exit -1
}

if [ $# -ne 1 ]; then usage; fi

SEQTEST=../../libs/seq_file/benchmarks/seqtest

# Load into disk cache
$SEQTEST --no-zlib --no-buf $1 >& /dev/null

(
time $SEQTEST --no-zlib --no-buf $1;
time $SEQTEST --no-zlib --no-buf $1;
time $SEQTEST --no-zlib --no-buf $1;
time $SEQTEST --no-zlib --no-buf $1;
time $SEQTEST --no-zlib --no-buf $1;

time $SEQTEST --no-zlib          $1;
time $SEQTEST --no-zlib          $1;
time $SEQTEST --no-zlib          $1;
time $SEQTEST --no-zlib          $1;
time $SEQTEST --no-zlib          $1;
) 2>&1 | \
grep '^user' | \
sed -E 's/.*([0-9]+)m([0-9\.]+)s.*/\1 \2/g' | \
awk '{print $1*60+$2}'
