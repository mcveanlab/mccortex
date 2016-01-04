#!/bin/bash

set -eou pipefail
set -o xtrace

CTXDIR=~/mccortex

VCF=/data2/users/turner/giab/NISTv2.19/GIAB2.19.vcf.gz
REF=/data1/users/turner/cortex_sims/NA12878/ref/human_g1k_v37.fasta.gz

REFGRAPH=../pipeline/k31/ref/ref.ctx
CLEANGRAPH=../pipeline/k31/graphs/NA12878.clean.ctx

VCFDIST=$CTXDIR/libs/vcf-slim/bin/vcfdist
VCFCONTIGS=$CTXDIR/libs/vcf-slim/bin/vcfcontigs
MCCORTEX=$CTXDIR/bin/mccortex31
CALL_COV_ANALYSIS=$CTXDIR/scripts/analysis/calls-covg.pl

$VCFDIST 31 $VCF > in.k31.vcf
$VCFCONTIGS --trim 31 $REF in.k31.vcf > in.k31.fa
$MCCORTEX coverage -m 50G --edges --degree --seq in.k31.fa --out in.k31.fa.covg $REFGRAPH $CLEANGRAPH >& in.k31.fa.covg.log
$CALL_COV_ANALYSIS in.k31.fa.covg >& coverage.txt
