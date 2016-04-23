#!/bin/bash

set -eou pipefail
set -o xtrace

#
# Calls were validated by mapping to PacBio assembly.
# Take validated calls and decompose with vcfallelicprimitives
# Then get SNP / indel counts.
# Create table: snp_indel.csv
# Create plot: kleb_indels.pdf
#

CTXDIR=~/mccortex
VCFHIST=$CTXDIR/libs/biogrok/vcf-hist
VCFSORT=$CTXDIR/libs/biogrok/vcf-sort
BGZIP=$CTXDIR/libs/htslib/bgzip
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools
VCFLIB_PRIMITIVES=~/Applications/homebrew/Cellar/vcflib/1.0.0/bin/vcfallelicprimitives

MCTXVCFS=../mcrun/mapping_truth/
PLATYPUSVCF=../platypus/mapping_truth/platypus.vcf.gz
FREEBAYESVCF=../freebayes/mapping_truth/freebayes.vcf.gz
CORTEXVCF=../cortex/mapping_truth/cortex.k31.k61.vcf.gz
BREAKPOINTSVCF=$MCTXVCFS/breakpoints.joint.links.k31.k61.vcf.gz
BUBBLESVCF=$MCTXVCFS/bubbles.joint.links.k31.k61.vcf.gz

mkdir -p vcfs

function decomp_vcf {
  local invcf=$1
  local outvcf=$(echo `dirname $2`/`basename $2 .gz`) # remove .gz
  [ -e $invcf.csi ] || $BCFTOOLS index $invcf
  $VCFLIB_PRIMITIVES $invcf | $VCFSORT > $outvcf
  $BGZIP $outvcf
  $BCFTOOLS index $outvcf.gz
}

decomp_vcf $PLATYPUSVCF vcfs/platypus.vcf.gz
decomp_vcf $FREEBAYESVCF vcfs/freebayes.vcf.gz
decomp_vcf $CORTEXVCF vcfs/cortex.vcf.gz
decomp_vcf $BREAKPOINTSVCF vcfs/breakpoints.vcf.gz
decomp_vcf $BUBBLESVCF vcfs/bubbles.vcf.gz

mkdir -p hists
for name in platypus freebayes cortex bubbles breakpoints; do
  $VCFHIST vcfs/$name.vcf.gz | tr ' ' ',' > hists/$name.hist.csv
done

./kleb-plot-hist.R

function snp_indel_counts {
  local csv=$1
  SNPS=$(grep '^0,' $csv | tr ',' '\t' | awk '{print $2}')
  INDELS=$(grep -v '^0,' $csv | awk -F, '{x = x+$2} END{print x}')
  printf "$csv\t$SNPS\t$INDELS\n"
}

# Get indel counts
(snp_indel_counts hists/platypus.hist.csv;
 snp_indel_counts hists/freebayes.hist.csv;
 snp_indel_counts hists/cortex.hist.csv;
 snp_indel_counts hists/bubbles.hist.csv;
 snp_indel_counts hists/breakpoints.hist.csv) > snp_indel.csv
