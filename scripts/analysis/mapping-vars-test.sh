#!/bin/bash

set -eou pipefail

if [ $# -ne 4 ]; then
  echo "Usage: $0 <in.vcf.gz> <ref.fa> <truth.fa> <out-dir>" 1>&2
  echo "  writes: OUT.fa, OUT.sam, OUT.stats.txt, OUT.sites.txt, OUT.sites.vcf.gz" 1>&2
  echo "  sites that map + pass are in: OUT.vcf.gz" 1>&2
  exit -1
fi
set -o xtrace

CTXDIR=$( cd $( dirname ${BASH_SOURCE[0]} ) && cd ../.. && pwd )
VCFCONTIGS=$CTXDIR/libs/vcf-slim/bin/vcfcontigs
SAM2VCF=$CTXDIR/libs/vcf-slim/scripts/sam-name-to-vcf.sh
BWA=$CTXDIR/libs/bwa/bwa
BGZIP=$CTXDIR/libs/htslib/bgzip
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools
SAMCMP=$CTXDIR/scripts/analysis/haploid-sam-compare.py

INVCF=$1
REF=$2
TRUTHFA=$3
PREFIX=$4

OUTFASTA=$PREFIX.fa
OUTSAM=$PREFIX.sam
OUTSTATS=$PREFIX.stats.txt
OUTSITES=$PREFIX.sites.txt
TMPVCF=$PREFIX.sites.vcf.gz
OUTVCF=$PREFIX.vcf.gz

mkdir -p $(dirname $OUTFASTA)

$VCFCONTIGS --max-alt 50 --trim --no-ref 50 $REF $INVCF > $OUTFASTA
$BWA mem $TRUTHFA $OUTFASTA > $OUTSAM
$SAMCMP --print-valid $OUTSITES $OUTSAM > $OUTSTATS
$SAM2VCF $OUTSITES | $BGZIP -c > $TMPVCF
$BCFTOOLS index $TMPVCF

$BCFTOOLS isec -n 2 -w 1 -o $OUTVCF -O z $INVCF $TMPVCF
$BCFTOOLS index $OUTVCF
