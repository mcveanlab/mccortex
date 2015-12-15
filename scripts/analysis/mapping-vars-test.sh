#!/bin/bash

set -eou pipefail

if [ $# -ne 5 ]; then
  echo "Usage: $0 <in.vcf> <ref.fa> <truth.fa> <out-fa> <out-sam>" 1>&2
  exit -1
fi
set -o xtrace

CTXDIR=$( cd $( dirname ${BASH_SOURCE[0]} ) && cd ../.. && pwd )
VCFCONTIGS=$CTXDIR/libs/vcf-slim/bin/vcfcontigs
BWA=$CTXDIR/libs/bwa/bwa
SAMCMP=$CTXDIR/scripts/analysis/haploid-sam-compare.py

INVCF=$1
REF=$2
TRUTHFA=$3
OUTFASTA=$4
OUTSAM=$5

$VCFCONTIGS --max-alt 50 --trim --no-ref 50 $REF $INVCF > $OUTFASTA
$BWA mem $TRUTHFA $OUTFASTA > $OUTSAM
$SAMCMP $OUTSAM
