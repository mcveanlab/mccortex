#!/bin/bash

CORTEXDIR=~/cortex/releases/CORTEX_release_v1.0.5.21
CTXDIR=../../../

RUNCALLS=$CORTEXDIR/scripts/calling/run_calls.pl
CORTEX=$CORTEXDIR/bin/cortex_var_31_c1
MCCORTEX=$CTXDIR/bin/mccortex31
STAMPY=/apps/well/stampy/1.0.23-py2.6/stampy.py
REF=$(readlink -f ../../data/chr22/chr22_17M_18M.fa)
VCFTOOLSDIR=~/bioinf/vcftools_0.1.12b/


$RUNCALLS \
--first_kmer 31 \
--last_kmer 31 \
--kmer_step 2 \
--fastaq_index samples.txt \
--auto_cleaning yes \
--bc yes \
--pd no \
--outdir cortex_run \
--outvcf chr22_17M_18M \
--ploidy 2 \
--stampy_hash ref/chr22_17M_18M \
--stampy_bin $STAMPY \
--list_ref_fasta ref/ref.falist \
--refbindir ref/ \
--genome_size 1000000 \
--qthresh 5 \
--mem_height 20 --mem_width 100 \
--vcftools_dir $VCFTOOLSDIR \
--do_union yes \
--ref CoordinatesAndInCalling \
--workflow independent \
--logfile runcalls.log
