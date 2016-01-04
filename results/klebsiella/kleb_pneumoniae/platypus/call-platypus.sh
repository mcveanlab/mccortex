#!/bin/bash

set -eou pipefail
set -o xtrace

REF=../ref/GCF_000016305.1_ASM1630v1_genomic.fa
BAM=../remap/mapped/KlebPneu.bam
BAMRMDUP=../remap/mapped/KlebPneu.rmdup.bam

CTXDIR=~/mccortex
VCF_ADD_CONTIGS=$CTXDIR/libs/biogrok/vcf-add-contigs
VCF_PASS=$CTXDIR/libs/biogrok/vcf-pass

PLATDIR=~/bioinf/Platypus

source $PLATDIR/prepare.sh
python $PLATDIR/bin/Platypus.py callVariants --logFileName platypus.rmdup.log \
                                             --output=platypus.rmdup.vcf \
                                             --refFile=$REF --bamFiles=$BAMRMDUP >& platypus.rmdup.vcf.log

python $PLATDIR/bin/Platypus.py callVariants --logFileName platypus.assem.log \
                                             --output=platypus.assem.vcf --assemble=1 \
                                             --refFile=$REF --bamFiles=$BAMRMDUP >& platypus.assem.vcf.log

# Add contigs to header
$VCF_ADD_CONTIGS <(dnacat --lengths $REF) KlebPneu_MGH_78578 platypus.rmdup.vcf | \
  $BCFTOOLS norm --check-ref x -m -any --fasta-ref $REF --site-win 5000 | \
  $BCFTOOLS norm --rm-dup any --do-not-normalize | \
  $VCF_PASS > platypus.vcf
$BGZIP platypus.vcf
$BCFTOOLS index platypus.vcf.gz

# Analysis
rm -rf mummer_isec mapping_truth cortex.k31.k61.{mapping,isec}.log
./analysis.sh >& analysis.log
