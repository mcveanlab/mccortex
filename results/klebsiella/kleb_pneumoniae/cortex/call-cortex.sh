#!/bin/bash

set -eou pipefail
set -o xtrace

#
# Cortex call at k=31,61
#

CWD=`pwd`
CTXDIR=~/mccortex
CORTEXDIR=~/bioinf/cortex
VCFTOOLS_DIR=~/bioinf/vcftools_0.1.12b
STAMPY_BIN=~/bioinf/stampy-1.0.23/stampy.py

VCF_ADD_CONTIGS=$CTXDIR/libs/biogrok/vcf-add-contigs
VCF_PASS=$CTXDIR/libs/biogrok/vcf-pass
TXT_SKIP_EMPTY=$CTXDIR/libs/biogrok/txt-skip-empty
BGZIP=$CTXDIR/libs/htslib/bgzip
BCFTOOLS=$CTXDIR/libs/bcftools/bcftools

function myreadlink() {
  ( cd $(dirname $1); echo $PWD/$(basename $1); )
}

REF=$(myreadlink ../ref/GCF_000016305.1_ASM1630v1_genomic.fa)
BAM=$(myreadlink ../remap/mapped/KlebPneu.bam)
SAMPLENAME=KlebPneu

# Set up environment
source $CORTEXDIR/prerun.sh

# Build ref
mkdir -p ref
echo $REF > ref/ref.falist

# ref k=31 (5.2GB)
time $CORTEXDIR/bin/cortex_var_31_c1 \
  --kmer 31 --se_list ref/ref.falist \
  --sample_id ref \
  --dump_binary ref/ref.k31.ctx \
  --mem_height 24 --mem_width 24 \
  >& ref/ref.k31.ctx.log

# ref k=61 (8.2GB)
time $CORTEXDIR/bin/cortex_var_63_c1 \
  --kmer 61 --se_list ref/ref.falist \
  --sample_id ref \
  --dump_binary ref/ref.k61.ctx \
  --mem_height 24 --mem_width 24 \
  >& ref/ref.k61.ctx.log

# Make stampy hash
cd ref
$STAMPY_BIN -G ref $REF
$STAMPY_BIN -g ref -H ref
cd ..

echo $BAM > sample.falist
printf "$SAMPLENAME\tsample.falist\t.\t.\n" > sample.index

# Run calling
# Estimate memory:
# ~/bioinf-perl/cortex_scripts/cortex_memory.pl -c 2 -k 31 24 24 => 7.1G
# ~/bioinf-perl/cortex_scripts/cortex_memory.pl -c 2 -k 63 24 24 => 10.1
# Actual memory (from top):
# k=31 => 
# k=61 => 
time $CORTEXDIR/scripts/calling/run_calls.pl \
  --auto_cleaning yes \
  --bc yes \
  --do_union yes \
  --first_kmer 31 \
  --last_kmer 61 \
  --kmer_step 30 \
  --genome_size 5694894 \
  --fastaq_index sample.index \
  --list_ref_fasta ref/ref.falist \
  --mem_height 24 \
  --mem_width 24 \
  --outdir cortex_run \
  --outvcf $SAMPLENAME.calls \
  --ploidy 1 \
  --qthresh 10 \
  --ref CoordinatesAndInCalling \
  --refbindir ref \
  --stampy_bin $STAMPY_BIN \
  --stampy_hash $CWD/ref/ref \
  --vcftools_dir $VCFTOOLS_DIR \
  --workflow joint \
  --logfile run.1.log


# Tidy up VCF
resultvcf=cortex_run/vcfs/KlebPneu.calls_wk_flow_J_RefCC_FINALcombined_BC_calls_at_all_k.decomp.vcf
$VCF_ADD_CONTIGS <(dnacat -L $REF) KlebPneu_MGH_78578 $resultvcf | \
  $TXT_SKIP_EMPTY | $VCF_PASS | \
  $BCFTOOLS norm --check-ref x -m -any --fasta-ref $REF --site-win 5000 | \
  $BCFTOOLS norm --rm-dup any --do-not-normalize > cortex.k31.k61.vcf
$BGZIP cortex.k31.k61.vcf
$BCFTOOLS index cortex.k31.k61.vcf.gz

# Analysis

rm -rf mummer_isec mapping_truth cortex.k31.k61.{mapping,isec}.log
./analysis.sh >& analysis.log

