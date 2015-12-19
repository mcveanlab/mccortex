#!/bin/bash

set -eou pipefail

if [ $# -ne 3 ]; then
  echo "Usage: $0 <mummer.vcf> <test.vcf> <isec-dir>" 1>&2
  exit -1
fi
set -o xtrace

TRUTHVCF=$1
TESTVCF=$2
ISECDIR=$3

CTXDIR=$( cd $( dirname ${BASH_SOURCE[0]} ) && cd ../.. && pwd )
VCFNLINES=$CTXDIR/libs/biogrok/vcf-count
VCFSORT=$CTXDIR/libs/biogrok/vcf-sort
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools
BGZIP=$CTXDIR/libs/htslib/bgzip
VCFLIB_PRIMITIVES=~/Applications/homebrew/Cellar/vcflib/1.0.0/bin/vcfallelicprimitives

function GIABisecstats {
  local isecdir=$1
  local TRUTH=`$VCFNLINES $isecdir/0000.vcf`
  local MCTX=`$VCFNLINES $isecdir/0001.vcf`
  local BOTH=`$VCFNLINES $isecdir/0002.vcf`
  local NTRUTH=$(($TRUTH+$BOTH))
  echo $TRUTH
  echo $MCTX
  echo $BOTH
  echo POW=`bc <<< "scale=3; $BOTH / $NTRUTH"`
  echo FP=`bc <<< "scale=3; $MCTX / ($BOTH+$MCTX)"`
}

function mk_isec {
  local TRUTH=$1
  local INVCF=$2
  local DIR=$3
  local TMPVCF=$DIR.vcf.gz
  if [ ! -e $TMPVCF ]; then
    # Break original VCF down into allelic primitives
    # take only SNPs
    [ -e $INVCF.csi ] || $BCFTOOLS index $INVCF
    $VCFLIB_PRIMITIVES $INVCF | $VCFSORT | \
      $BCFTOOLS view --types snps --output-file $TMPVCF --output-type z -
    $BCFTOOLS index $TMPVCF
  fi
  if [ ! -d $DIR ]; then
    mkdir -p $DIR
    $BCFTOOLS isec $TRUTH $TMPVCF -p $DIR
  fi
}

mk_isec $TRUTHVCF $TESTVCF $ISECDIR
GIABisecstats $ISECDIR
