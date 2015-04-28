#!/bin/bash

set -euo pipefail

CORTEXDIR=~/cortex/releases/CORTEX_release_v1.0.5.21
CTXDIR=../../../

RUNCALLS=$CORTEXDIR/scripts/calling/run_calls.pl
CORTEX=$CORTEXDIR/bin/cortex_var_31_c1
MCCORTEX=$CTXDIR/bin/mccortex31
STAMPY=/apps/well/stampy/1.0.23-py2.6/stampy.py
REF=$(readlink -f ../../data/chr22/chr22_17M_18M.fa)
VCFTOOLSDIR=~/bioinf/vcftools_0.1.12b/

mkdir -p ref
echo $REF > ref/ref.falist

# Make stampy hash
$STAMPY -G ref/chr22_17M_18M $REF
$STAMPY -g ref/chr22_17M_18M -H ref/chr22_17M_18M

# Build reference graph file
#$CORTEX --mem_height 20 --mem_width 100 --kmer_size 31 --sample_name REF --se_list ref/ref.falist --dump_binary ref/ref.k31.ctx
$MCCORTEX build -k 31 -s REF -1 $REF ref/ref.k31.ctx >& ref/ref.k31.ctx.log

(readlink -f ../reads/chrom0.30X.1.fa.gz;
 readlink -f ../reads/chrom1.30X.1.fa.gz) > reads.1.falist
(readlink -f ../reads/chrom0.30X.2.fa.gz;
 readlink -f ../reads/chrom1.30X.2.fa.gz) > reads.2.falist

printf "MrSample\t.\t%s\t%s\n" reads.1.falist reads.2.falist > samples.txt
