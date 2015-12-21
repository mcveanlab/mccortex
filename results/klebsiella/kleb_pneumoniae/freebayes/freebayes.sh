#!/bin/bash

set -eou pipefail
set -o xtrace

REF=/data2/users/turner/cortex_sims/klebsiella/kleb_pneumoniae/ref/GCF_000016305.1_ASM1630v1_genomic.fa
BAM=/data2/users/turner/cortex_sims/klebsiella/kleb_pneumoniae/remap/mapped/KlebPneu.bam
BAMRMDUP=/data2/users/turner/cortex_sims/klebsiella/kleb_pneumoniae/remap/mapped/KlebPneu.rmdup.bam

CTXDIR=~/mccortex
FREEBAYES=~/bioinf/freebayes/bin/freebayes
BGZIP=$CTXDIR/libs/htslib/bgzip
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools

$FREEBAYES -f $REF -p 1 $BAMRMDUP > freebayes.rmdup.vcf
$BGZIP freebayes.rmdup.vcf
$BCFTOOLS index freebayes.rmdup.vcf.gz

$BCFTOOLS norm --check-ref x -m -any --fasta-ref $REF --site-win 5000 freebayes.rmdup.vcf.gz | \
  $BCFTOOLS norm --rm-dup any --do-not-normalize | \
  $VCF_PASS > freebayes.vcf
$BGZIP freebayes.vcf
$BCFTOOLS index freebayes.vcf.gz

# Analysis
rm -rf mummer_isec mapping_truth cortex.k31.k61.{mapping,isec}.log
./analysis.sh >& analysis.log
