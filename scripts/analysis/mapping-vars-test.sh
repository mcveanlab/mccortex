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
BWA=$CTXDIR/libs/bwa/bwa
BGZIP=$CTXDIR/libs/htslib/bgzip
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools
VCFCONTIGS=$CTXDIR/libs/vcf-slim/bin/vcfcontigs
SAM2VCF=$CTXDIR/libs/vcf-slim/scripts/sam-name-to-vcf.sh
VCFRENAME=$CTXDIR/libs/biogrok/vcf-rename
VCF_SELECT_ID=$CTXDIR/libs/biogrok/vcf-select-id
SAMCMP=$CTXDIR/scripts/analysis/haploid-sam-compare.py


INVCF=$1
REF=$2
TRUTHFA=$3
PREFIX=$4

OUTFASTA=$PREFIX.fa
OUTSAM=$PREFIX.sam
OUTSTATS=$PREFIX.stats.txt
OUTSITES=$PREFIX.sites.txt
RENAMEDVCF=$PREFIX.renamed.vcf.gz
OUTVCF=$PREFIX.vcf.gz

mkdir -p $(dirname $OUTFASTA)

$VCFRENAME $INVCF > $RENAMEDVCF
$VCFCONTIGS --trim --no-ref 50 $REF $RENAMEDVCF > $OUTFASTA
$BWA mem $TRUTHFA $OUTFASTA > $OUTSAM
$SAMCMP --print-valid $OUTSITES $OUTSAM > $OUTSTATS
$VCF_SELECT_ID <(cut -d: -f3 $OUTSITES) $RENAMEDVCF | $BGZIP -c > $OUTVCF
