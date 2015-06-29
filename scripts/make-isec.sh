#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail

if [[ $# -ne 3 ]]
then
  echo "usage: $0 <tmpdir> <truth.vcf.gz> <results.vcf.gz>" 1>&2
  echo "  Create tmp dir and use it to count intersection of indexed VCFs" 1>&2
  exit -1
fi

TMPDIR="$1"
TRUTHVCF="$2"
RESULTSVCF="$3"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CTXDIR="$DIR/.."

BCFTOOLS=$CTXDIR"/libs/bcftools/bcftools"
VCFALLELES=$CTXDIR"/libs/biogrok/vcf-count-alleles"

$BCFTOOLS isec $TRUTHVCF $RESULTSVCF -p $TMPDIR

MISSED=`$VCFALLELES $TMPDIR/0000.vcf`
FP=`$VCFALLELES $TMPDIR/0001.vcf`
FOUND=`$VCFALLELES $TMPDIR/0002.vcf`
NCALLED=`$VCFALLELES $RESULTSVCF`
NTRUTH=`$VCFALLELES $TRUTHVCF`

awk 'BEGIN{printf("Missed: %4d / %4d (%5.2f%%)\n",'$MISSED','$NTRUTH',100*'$MISSED'/'$NTRUTH')}'
awk 'BEGIN{printf("FP:     %4d / %4d (%5.2f%%)\n",'$FP','$NCALLED',100*'$FP'/'$NCALLED')}'
awk 'BEGIN{printf("Found:  %4d / %4d (%5.2f%%)\n",'$FOUND','$NTRUTH',100*'$FOUND'/'$NTRUTH')}'

echo "remember to delete temp dir: $TMPDIR" 1>&2
