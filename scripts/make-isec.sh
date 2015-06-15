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
A=`$VCFALLELES $TMPDIR/0000.vcf`
B=`$VCFALLELES $TMPDIR/0001.vcf`
C=`$VCFALLELES $TMPDIR/0002.vcf`
X=`$VCFALLELES $RESULTSVCF`
T=`$VCFALLELES $TRUTHVCF`
awk 'BEGIN{printf("Missed: %4d / %4d (%5.2f%%)\n",'$A','$T',100*'$A'/'$T')}'
awk 'BEGIN{printf("FP:     %4d / %4d (%5.2f%%)\n",'$B','$X',100*'$B'/'$X')}'
awk 'BEGIN{printf("Found:  %4d / %4d (%5.2f%%)\n",'$C','$T',100*'$C'/'$T')}'

echo "remember to delete temp dir: $TMPDIR" 1>&2
