#!/bin/bash

set -euo pipefail

CTXDIR=../../..
VCFHEADER=$CTXDIR/scripts/bash/vcf-header
#VCFSORT=$CTXDIR/scripts/bash/vcf-sort
VCFSORT=./vcf-sort
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools
BGZIP=$CTXDIR/libs/htslib/bgzip
VCFREF=~/c/vcf-hack/bin/vcfref

# Files
REF=../../data/chr22/chr22_17M_18M.fa
IN=cortex_run/vcfs/chr22_17M_18M_union_BC_calls_k31.decomp.vcf

# Add '##contig=<ID=chr22_17M_18M,length=1000000,assembly=hg19>'
# to header, and fix an INFO field
( $VCFHEADER $IN | \
  grep -v '^##contig' | \
  grep -v '^#CHROM' | \
  sed 's/, Description=/,Description=/g';
  echo '##INFO=<ID=KMER,Number=1,Type=Integer,Description="Kmer used for calling">';
  echo '##contig=<ID=chr22_17M_18M,length=1000000,assembly=hg19>';
  $VCFHEADER $IN | grep '^#CHROM' ) > new_header.txt

# Put new header on and filter ref mismatches
( cat new_header.txt;
  $VCFREF -s $IN $REF | grep -v '^#' |  sort -k1,1d -k2,2n ) > cortex.sort.vcf

$BCFTOOLS norm --remove-duplicates --fasta-ref $REF --multiallelics +both cortex.sort.vcf > cortex.norm.vcf
$BGZIP cortex.norm.vcf
$BCFTOOLS index cortex.norm.vcf.gz
