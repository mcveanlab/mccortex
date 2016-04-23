#!/bin/bash

set -eou pipefail
set -o xtrace

CTXDIR=../../../../

function myreadlink() {
  ( cd $(dirname $1); echo $PWD/$(basename $1); )
}

REF=$(myreadlink ../ref/GCF_000016305.1_ASM1630v1_genomic.fna.gz)
MUMMER=$(myreadlink ../mummer/mummer.vcf.gz)
TRUTH=$(myreadlink ../truth/CAV1016.fa)
BCFTOOLS=$(myreadlink $CTXDIR/libs/bcftools/bcftools)
MAPPING_TEST=$(myreadlink $CTXDIR/scripts/analysis/mapping-vars-test.sh)
MUMMER_ISEC=$(myreadlink $CTXDIR/scripts/analysis/mummer-vcf-isec.sh)

mkdir -p mapping_truth mummer_isec

for vcf in `ls vcfs/*k{61,51}.vcf.gz`; do
  name=`basename $vcf .vcf.gz`
  echo "== $name"
  [ -e $vcf.csi ] || $BCFTOOLS index $vcf
  $MAPPING_TEST $vcf $REF $TRUTH mapping_truth/$name
  $MUMMER_ISEC $MUMMER $vcf mummer_isec/$name >& $name.isec.log
done

