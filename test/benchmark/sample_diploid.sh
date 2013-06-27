#!/bin/bash

# Set exit on error
set -e

SCRIPT=$0

if [[ $# -ne 8 ]]
then
	echo "Usage: $0 <genome.fa> <snps> <indels> <invs> <invlen> <readlen> <mpsize> <covg>"
	exit 1
fi

INPUT_SEQ=$1
SNPS=$2
INDELS=$3
INV=$4
INVLEN=$5
READLEN=$6
MPSIZE=$7
ALLELECOVG=$8

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CTX="~/cortex/versions/current/bin/cortex_var_31_c1_s8 --kmer_size 31 --mem_width 8 --mem_height 18"
THREAD="$DIR/../../bin/ctx_thread_k31_c1"
CTX2="$DIR/../../bin/ctx_call_k31_c1"
PROC="$DIR/../../bin/ctx_unique"
PLACE="$DIR/../../bin/ctx_place"
BIOINF="$DIR/../../libs/bioinf-perl"
# BIOINF=~/bioinf-perl
HAPLEN="$DIR/longest-haplotype.sh"

if [[ ! -e $INPUT_SEQ ]]
then
	echo "$0: cannot read $INPUT_SEQ"; exit 1;
fi

cmd()
{
	echo $@
	eval "$@" || { echo "$SCRIPT: Command failed: '$@'"; exit 1; }
}

# Simulate mutations
cmd "zcat -f $INPUT_SEQ | $BIOINF/sim_mutations/sim_mutations.pl --snps $SNPS --indels $INDELS --invs $INV --invlen $INVLEN ./ 2 -"

# Generate reads
cmd $BIOINF/sim_mutations/sim_reads.pl --readlen $READLEN --mpsize $MPSIZE --covg $ALLELECOVG reads0 genome0.fa
cmd $BIOINF/sim_mutations/sim_reads.pl --readlen $READLEN --mpsize $MPSIZE --covg $ALLELECOVG reads1 genome1.fa

cmd "echo reads0.0.fa > reads.0.falist"
cmd "echo reads1.0.fa >> reads.0.falist"

cmd "echo reads0.1.fa > reads.1.falist"
cmd "echo reads1.1.fa >> reads.1.falist"

cmd "echo reads0.0.fa > reads.se.falist"
cmd "echo reads0.1.fa >> reads.se.falist"
cmd "echo reads1.0.fa >> reads.se.falist"
cmd "echo reads1.1.fa >> reads.se.falist"

# Make diploid.k31.ctx
cmd $CTX --pe_list reads.0.falist,reads.1.falist --dump_binary diploid.k31.ctx

# Add shades to binary
# cmd $CTX --load_binary diploid.k31.ctx --add_shades --pe_list reads.0.falist,reads.1.falist --dump_binary diploid.shaded.k31.ctx

# Call with old bc
cmd time $CTX --load_binary diploid.k31.ctx --detect_bubbles1 0/0 --output_bubbles1 diploid.oldbc.bubbles --print_colour_coverages
cmd $PROC diploid.oldbc.bubbles diploid.oldbc
cmd gzip -d -f diploid.oldbc.vcf.gz

# Fix buggy output from old bc
grep -v 'LF=;' diploid.oldbc.vcf > diploid.oldbc.vcf.2;
mv diploid.oldbc.vcf.2 diploid.oldbc.vcf

# Call with new bc
cmd time $CTX --load_binary diploid.k31.ctx --paths_caller diploid.newbc.bubbles.gz
cmd $PROC diploid.newbc.bubbles.gz diploid.newbc
cmd gzip -d -f diploid.newbc.vcf.gz

# Call with new bc + shades (also add shades)
# cmd time $CTX --load_binary diploid.shaded.k31.ctx --paths_caller diploid.newbc.shaded.bubbles.gz
cmd time $CTX --load_binary diploid.k31.ctx --add_shades --pe_list reads.0.falist,reads.1.falist --paths_caller diploid.newbc.shaded.bubbles.gz
cmd $PROC diploid.newbc.shaded.bubbles.gz diploid.newbc.shaded
cmd gzip -d -f diploid.newbc.shaded.vcf.gz

# Call with ctx2
# cmd time $THREAD --pe_list 0 reads.0.falist reads.1.falist 2 100MB diploid.k31.ctx
cmd time $THREAD --se_list 0 reads.se.falist 2 100M diploid.k31.ctx
# cmd time $THREAD --se_list 0 reads.se.falist --pe_list 0 reads.0.falist reads.1.falist 100M diploid.k31.ctx
cmd time $CTX2 diploid.k31.ctx 100M diploid.pac.bubbles.gz
cmd $PROC diploid.pac.bubbles.gz diploid.pac
cmd gzip -d -f diploid.pac.vcf.gz

# Generate truth VCF
cmd "$BIOINF/sim_mutations/sim_vcf.pl 31 genome0.fa mask0.fa genome1.fa mask1.fa > truth.k31.vcf"

# Compare
cmd $BIOINF/sim_mutations/sim_compare.pl truth.k31.vcf diploid.oldbc.vcf truth.k31.oldbc.vcf OLDBC falsepos.k31.oldbc.vcf genome0.fa genome1.fa
cmd $HAPLEN diploid.oldbc.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl truth.k31.oldbc.vcf diploid.newbc.vcf truth.k31.newbc.vcf NEWBC falsepos.k31.newbc.vcf genome0.fa genome1.fa
cmd $HAPLEN diploid.newbc.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl truth.k31.newbc.vcf diploid.newbc.shaded.vcf truth.k31.shaded.vcf SHADED falsepos.k31.shaded.vcf genome0.fa genome1.fa
cmd $HAPLEN diploid.newbc.shaded.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl truth.k31.shaded.vcf diploid.pac.vcf truth.k31.pac.vcf PAC falsepos.k31.pac.vcf genome0.fa genome1.fa
cmd $HAPLEN diploid.pac.vcf
cmd $HAPLEN truth.k31.vcf

exit

# Map 5' flanks
STAMPY=stampy.py
if [ ! -e stampy.py ]; then STAMPY="python2.6 $HOME/bioinf/stampy-1.0.20/stampy.py"; fi
if [ ! -e ../chr21.stidx ]; then `$STAMPY -G ../chr21 ../chr21.1Mb.fa.gz`; fi
if [ ! -e ../chr21.sthash ]; then `$STAMPY -g ../chr21 -H ../chr21`; fi

cmd "$STAMPY -g ../chr21 -h ../chr21 --inputformat=fasta -M diploid.oldbc.5pflanks.fa.gz > diploid.oldbc.5pflanks.sam"
cmd "$STAMPY -g ../chr21 -h ../chr21 --inputformat=fasta -M diploid.newbc.5pflanks.fa.gz > diploid.newbc.5pflanks.sam"
cmd "$STAMPY -g ../chr21 -h ../chr21 --inputformat=fasta -M diploid.newbc.shaded.5pflanks.fa.gz > diploid.newbc.shaded.5pflanks.sam"
cmd "$STAMPY -g ../chr21 -h ../chr21 --inputformat=fasta -M diploid.pac.5pflanks.fa.gz > diploid.pac.5pflanks.sam"

# Place calls
cmd "time $PLACE diploid.oldbc.vcf diploid.oldbc.5pflanks.sam ../chr21.1Mb.fa.gz > diploid.oldbc.decomp.vcf"
cmd "time $PLACE diploid.newbc.vcf diploid.newbc.5pflanks.sam ../chr21.1Mb.fa.gz > diploid.newbc.decomp.vcf"
cmd "time $PLACE diploid.newbc.shaded.vcf diploid.newbc.shaded.5pflanks.sam ../chr21.1Mb.fa.gz > diploid.newbc.shaded.decomp.vcf"
cmd "time $PLACE diploid.pac.vcf diploid.pac.5pflanks.sam ../chr21.1Mb.fa.gz > diploid.pac.decomp.vcf"
