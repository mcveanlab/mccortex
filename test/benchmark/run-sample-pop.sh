#!/bin/bash

# Set exit on error
set -e

SCRIPT=$0

if [[ $# -ne 11 ]]
then
	echo "Usage: $0 <genome.fa> <indivs> <ploidy> <kmer> <snps> <indels> <invs> <invlen> <readlen> <mpsize> <covg>"
	exit 1
fi

INPUT_SEQ=$1
NUM_INDIVS=$2
PLOIDY=$3
KMER=$4
SNPS=$5
INDELS=$6
INV=$7
INVLEN=$8
READLEN=$9
MPSIZE=${10}
ALLELECOVG=${11}

TLEN=$READLEN+$MPSIZE+$READLEN

NCHROMS=$(($NUM_INDIVS * $PLOIDY))

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OLDCTX="~/cortex/releases/CORTEX_release_v1.0.5.20/bin/cortex_var_31_c$NUM_INDIVS --kmer_size $KMER --mem_height 18 --mem_width 20"
BUILDCTX="$DIR/../../bin/ctx31 build"
THREAD="$DIR/../../bin/ctx31 thread"
CALL="$DIR/../../bin/ctx31 call"
PROC="$DIR/../../bin/ctx31 unique"
PLACE="$DIR/../../bin/ctx31 place"
BIOINF="$DIR/../../libs/bioinf-perl"
HAPLEN="$DIR/longest-haplotype.sh"
OLDCLEAN="$DIR/clean_bubbles.pl"
READSIM="$DIR../../libs/readsim"


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
cmd "zcat -f $INPUT_SEQ | $BIOINF/sim_mutations/sim_mutations.pl --snps $SNPS --indels $INDELS --invs $INV --invlen $INVLEN ./ $NCHROMS -"

LASTCHROM=$(($NCHROMS-1))
LASTINDIV=$(($NUM_INDIVS-1))

# Generate reads
for i in $(seq 0 $LASTCHROM)
do
  # cmd $READSIM -r genome$i.fa -t $TLEN -v 0 -l $READLEN -d $ALLELECOVG ../emptyprofile.fq reads$i
  cmd $BIOINF/sim_mutations/sim_reads.pl --readlen $READLEN --mpsize $MPSIZE --covg $ALLELECOVG reads$i genome$i.fa
  cmd "echo reads$i.0.fq.gz >> reads.0.falist"
  cmd "echo reads$i.1.fq.gz >> reads.1.falist"
done

# Make diploid.ctx
for i in $(seq 0 $LASTINDIV)
do
  a=$(($i*2))
  b=$(($i*2+1))
  PELIST="--seq2 reads$a.0.fa reads$a.1.fa --seq2 reads$b.0.fa reads$b.1.fa"
  cmd time $BUILDCTX -k $KMER -m 100MB --sample MrDiploid$i $PELIST diploid$i.ctx
done

# Merge
cmd ctx31 join -m 100M pop.ctx diploid{0..$LASTINDIV}.ctx

# Call with old bc
cmd time $OLDCTX --multicolour_bin pop.ctx --detect_bubbles1 -1/-1 --output_bubbles1 diploid.oldbc.bubbles --print_colour_coverages
# Fix buggy output from old bc
cmd mv diploid.oldbc.bubbles diploid.oldbc.bubbles.2
cmd "$OLDCLEAN $KMER diploid.oldbc.bubbles.2 > diploid.oldbc.bubbles"
cmd $PROC diploid.oldbc.bubbles diploid.oldbc
cmd gzip -d -f diploid.oldbc.vcf.gz

# Fix buggy output from old bc
# grep -v 'LF=;' diploid.oldbc.vcf > diploid.oldbc.vcf.2;
# mv diploid.oldbc.vcf.2 diploid.oldbc.vcf

# Call with new bc
cmd time $CALL -t 1 pop.ctx diploid.newbc.bubbles.gz
cmd $PROC diploid.newbc.bubbles.gz diploid.newbc
cmd gzip -d -f diploid.newbc.vcf.gz

# Call with new bc + shades (also add shades)
# cmd time $OLDCTX --load_binary pop.ctx --add_shades --pe_list reads.0.falist,reads.1.falist --paths_caller diploid.newbc.shaded.bubbles.gz
# cmd $PROC diploid.newbc.shaded.bubbles.gz diploid.newbc.shaded
# cmd gzip -d -f diploid.newbc.shaded.vcf.gz

SELIST=
PELIST=
SEPELIST=

for i in $(seq 0 $LASTINDIV)
do
  a=$(($i*2))
  b=$(($i*2+1))
  se="--seq reads$a.0.fa --seq reads$a.1.fa --seq reads$b.0.fa --seq reads$b.1.fa"
  pe="--seq2 reads$a.0.fa reads$a.1.fa --seq2 reads$b.0.fa reads$b.1.fa"
  SELIST="$SELIST --col $i $i $se"
  PELIST="$PELIST --col $i $i $pe"
  SEPELIST="$SEPELIST --col $i $i $se $pe"
done

cmd time $THREAD -t 1 $SELIST $NUM_INDIVS pop.se.ctp pop.ctx pop.ctx
cmd time $THREAD -t 1 $PELIST $NUM_INDIVS pop.pe.ctp pop.ctx pop.ctx
cmd time $THREAD -t 1 $SEPELIST $NUM_INDIVS pop.sepe.ctp pop.ctx pop.ctx

for x in se pe sepe
do
  cmd time $CALL -t 1 -m 100MB -p pop.$x.ctp pop.ctx diploid.pac.$x.bubbles.gz
  cmd $PROC diploid.pac.$x.bubbles.gz diploid.pac.$x
  cmd gzip -d -f diploid.pac.$x.vcf.gz
done

# Generate truth VCF
# MG is list of Mask+Genome files
MGLIST=
for i in $(seq 0 $LASTCHROM)
do
  MGLIST="$MGLIST genome$i.fa mask$i.fa"
done

cmd "$BIOINF/sim_mutations/sim_vcf.pl $KMER $MGLIST > truth.vcf"

# Compare
LVCF=truth.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF diploid.oldbc.vcf truth.oldbc.vcf OLDBC falsepos.oldbc.vcf genome{0..$LASTCHROM}.fa
cmd $HAPLEN diploid.oldbc.vcf
LVCF=truth.oldbc.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF diploid.newbc.vcf truth.newbc.vcf NEWBC falsepos.newbc.vcf genome{0..$LASTCHROM}.fa
cmd $HAPLEN diploid.newbc.vcf
LVCF=truth.newbc.vcf
# cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF diploid.newbc.shaded.vcf truth.shaded.vcf SHADED falsepos.shaded.vcf genome{0..$LASTCHROM}.fa
# cmd $HAPLEN diploid.newbc.shaded.vcf
# LVCF=truth.shaded.vcf

for x in se pe sepe
do
  cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF diploid.pac.$x.vcf truth.pac.$x.vcf PAC falsepos.pac.$x.vcf genome{0..$LASTCHROM}.fa
  cmd $HAPLEN diploid.pac.$x.vcf
  LVCF=truth.pac.$x.vcf
done

# cmd $HAPLEN truth.vcf

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
