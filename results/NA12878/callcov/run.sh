#!/bin/bash

set -eou pipefail
set -o xtrace

CTXDIR=~/mccortex

VCF=/data2/users/turner/giab/NISTv2.19/GIAB2.19.vcf.gz
REF=/data1/users/turner/cortex_sims/NA12878/ref/human_g1k_v37.fasta.gz

REFGRAPH=../pipeline/k31/ref/ref.ctx
CLEANGRAPH=../pipeline/k31/graphs/NA12878.clean.ctx

MCCORTEX=$CTXDIR/bin/mccortex31
VCFDIST=$CTXDIR/libs/vcf-slim/bin/vcfdist
VCFCONTIGS=$CTXDIR/libs/vcf-slim/bin/vcfcontigs
CALL_COV_ANALYSIS=$CTXDIR/scripts/analysis/calls-covg.pl
BGZIP=$CTXDIR/libs/htslib/bgzip
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools
VCFCOUNT=$CTXDIR/libs/biogrok/vcf-count

# $VCFDIST 31 $VCF > in.k31.vcf
# $VCFCONTIGS --trim 31 $REF in.k31.vcf > in.k31.fa
# $MCCORTEX coverage -m 100G --edges --degree --seq in.k31.fa --out in.k31.fa.covg $REFGRAPH $CLEANGRAPH >& in.k31.fa.covg.log
# $CALL_COV_ANALYSIS in.k31.fa.covg >& coverage.txt

# Find intersection to get power / FP for selected variants

function isec_stats {
  isecdir=$1
  TRTH=`$VCFCOUNT $isecdir/0000.vcf`
  MCTX=`$VCFCOUNT $isecdir/0001.vcf`
  BOTH=`$VCFCOUNT $isecdir/0002.vcf`
  TOTAL=$(($TRTH+$BOTH))
  echo $TRTH
  echo $MCTX
  echo $BOTH
  echo POW=`bc <<< "scale=3; $BOTH / $TOTAL"`
  echo FP=`bc <<< "scale=3; $MCTX / ($BOTH+$MCTX)"`
}

function get_stats {
  truthvcf=$1
  samplevcf=$2
  isecdir=$3
  echo "-- $truthvcf $samplevcf $isecdir/"
  if [[ ! -d $isecdir ]]; then
    bcftools isec $truthvcf $samplevcf -p $isecdir
  fi
  isec_stats $isecdir
}

#$BGZIP in.k31.vcf
#$BCFTOOLS index in.k31.vcf.gz

get_stats in.k31.vcf.gz ../rerun_decomp/k31/bubbles_plain/joint.bub.norm.vcf.gz isec_k31bub_plain
get_stats in.k31.vcf.gz ../rerun_decomp/k31/bubbles_links/joint.bub.norm.vcf.gz isec_k31bub_links
get_stats in.k31.vcf.gz ../rerun_decomp/k31/breakpoints_plain/joint.brk.norm.vcf.gz isec_k31brk_plain
get_stats in.k31.vcf.gz ../rerun_decomp/k31/breakpoints_links/joint.brk.norm.vcf.gz isec_k31brk_links

