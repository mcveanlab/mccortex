#!/bin/bash

# Set exit on error
set -e
# Use aliases
shopt -s expand_aliases

alias vcf-pass="awk '\$1 ~ /^#/ || \$7 ~ /^(PASS|\.)$/'"

cmd()
{
  echo $@
  if [ $RUN == "run" ]
  then
    eval "$@" || { echo "$SCRIPT: Command failed: '$@'"; exit 1; }
  fi
}

msg()
{
  if [ $RUN == "run" ]
  then
    echo $@
  else
    echo echo $@
  fi
}

SCRIPT=$0

if [[ $# -ne 14 ]]
then
	echo "Usage: $0 <genome.fa> <stampyhsh> <indivs> <ploidy> <kmer> <snps> 
                  <indels> <invs> <invlen> <readlen> <mpsize> <allelecovg>
                  <error|noerror> <run|print>"
	exit 1
fi

INPUT_SEQ=$1
STAMPY_HSH=$2
NUM_INDIVS=$3
PLOIDY=$4
KMER=$5
SNPS=$6
INDELS=$7
INV=$8
INVLEN=$9
READLEN=${10}
MPSIZE=${11}
ALLELECOVG=${12}
witherror=${13}
RUN=${14}

if [ $witherror != "error" ] && [ $witherror != "noerror" ]; then
  echo "Must specify <error|noerror>"; exit 1;
fi

if [ $RUN != "run" ] && [ $RUN != "print" ]; then
  echo "Must specify <run|print>"; exit 1;
fi

NCHROMS=$(($NUM_INDIVS * $PLOIDY))

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RELEASECTX="~/cortex/releases/CORTEX_release_v1.0.5.20/bin/cortex_var_31_c$NUM_INDIVS --kmer_size $KMER --mem_height 18 --mem_width 20"
SHADECTX="~/cortex/versions/current/bin/cortex_var_31_c"$NUM_INDIVS"_s8 --kmer_size $KMER --mem_height 18 --mem_width 20"
BUILDCTX="$DIR/../../bin/ctx31 build"
CLEANCTX="$DIR/../../bin/ctx31 clean"
JOINCTX="$DIR/../../bin/ctx31 join"
INFERCTX="$DIR/../../bin/ctx31 inferedges --pop"
THREADCTX="$DIR/../../bin/ctx31 thread"
CALLCTX="$DIR/../../bin/ctx31 call"
PROCCTX="$DIR/../../bin/ctx31 unique"
PLACECTX="$DIR/../../bin/ctx31 place"
TRAVERSE="$DIR/../../bin/traversal31"

BIOINF="$DIR/../../libs/bioinf-perl"
HAPLEN="$DIR/longest-haplotype.sh"
OLDCLEAN="$DIR/clean_bubbles.pl"
READSIM="$DIR/../../libs/readsim/readsim"
BCFTOOLS="~/bioinf/bcftools/bcftools"

CALIB="$DIR/PhiX.100K.1.fq.gz"


if [[ ! -e $INPUT_SEQ ]]
then
	echo "$0: warning cannot read $INPUT_SEQ" 1>&2
fi

if [ $RUN == "print" ]
then
  echo "#!/bin/bash"
  echo set -e
  echo shopt -s expand_aliases
  echo 'alias vcf-pass="awk '"'"'\$1 ~ /^#/ || \$7 ~ /^(PASS|\.)$/'"'"'"'
fi

cmd mkdir -p genomes reads graphs bubbles vcfs

# Simulate mutations
cmd "zcat -f $INPUT_SEQ | $BIOINF/sim_mutations/sim_mutations.pl --snps $SNPS --indels $INDELS --invs $INV --invlen $INVLEN genomes/ $NCHROMS -"

LASTCHROM=$(($NCHROMS-1))
LASTINDIV=$(($NUM_INDIVS-1))

if [ $witherror == "error" ]
then
  USECALIB="-p $CALIB"
fi

# Generate reads
for i in $(seq 0 $LASTCHROM)
do
  cmd "$READSIM -r <(cat genomes/genome$i.fa | tr -d '-') -i $MPSIZE -v 0.2 -l $READLEN -d $ALLELECOVG $USECALIB reads/reads$i"
  cmd "echo reads$i.1.fa.gz >> reads/reads$i.1.falist"
  cmd "echo reads$i.2.fa.gz >> reads/reads$i.2.falist"
done

# Make diploid.ctx
for i in $(seq 0 $LASTINDIV)
do
  a=$(($i*2))
  b=$(($i*2+1))
  PELIST="--seq2 reads/reads$a.1.fa.gz reads/reads$a.2.fa.gz --seq2 reads/reads$b.1.fa.gz reads/reads$b.2.fa.gz"
  cmd time $BUILDCTX -k $KMER -m 100MB --sample MrDiploid$i $PELIST graphs/diploid$i.ctx

  if [ $witherror == "error" ]
  then
    cmd time $CLEANCTX graphs/diploid$i.clean.ctx graphs/diploid$i.ctx
  fi
done

ctxext=ctx
if [ $witherror == "error" ]
then
  ctxext=clean.ctx
fi

# Merge
cmd time $JOINCTX -m 100M graphs/pop.ctx graphs/diploid{0..$LASTINDIV}.$ctxext
cmd time $INFERCTX graphs/pop.ctx

# Call with old bc
cmd time $RELEASECTX --multicolour_bin graphs/pop.ctx --detect_bubbles1 -1/-1 --output_bubbles1 bubbles/diploid.oldbc.bubbles --print_colour_coverages
# Fix buggy output from old bc
cmd mv bubbles/diploid.oldbc.bubbles bubbles/diploid.oldbc.bubbles.2
cmd "$OLDCLEAN $KMER bubbles/diploid.oldbc.bubbles.2 > bubbles/diploid.oldbc.bubbles"
cmd $PROCCTX bubbles/diploid.oldbc.bubbles vcfs/diploid.oldbc
cmd gzip -d -f vcfs/diploid.oldbc.vcf.gz

# Call with new bc
cmd time $CALLCTX -t 1 graphs/pop.ctx bubbles/diploid.newbc.bubbles.gz
cmd $PROCCTX bubbles/diploid.newbc.bubbles.gz vcfs/diploid.newbc
cmd gzip -d -f vcfs/diploid.newbc.vcf.gz

SELIST=
PELIST=
SEPELIST=
SHADEDLIST=

for i in $(seq 0 $LASTINDIV)
do
  a=$(($i*2))
  b=$(($i*2+1))
  se="--seq reads/reads$a.1.fa.gz --seq reads/reads$a.2.fa.gz --seq reads/reads$b.1.fa.gz --seq reads/reads$b.2.fa.gz"
  pe="--seq2 reads/reads$a.1.fa.gz reads/reads$a.2.fa.gz --seq2 reads/reads$b.1.fa.gz reads/reads$b.2.fa.gz"
  SELIST="$SELIST --col $i $i $se"
  PELIST="$PELIST --col $i $i $pe"
  SEPELIST="$SEPELIST --col $i $i $se $pe"
  SHADEDLIST="$SHADEDLIST --pe_list reads/reads$i.1.falist,reads/reads$i.2.falist"
done

# Call with new bc + shades (also add shades)
cmd time $SHADECTX --load_binary graphs/pop.ctx --add_shades $SHADEDLIST --paths_caller bubbles/diploid.shaded.bubbles.gz --paths_caller_cols -1
cmd $PROCCTX bubbles/diploid.shaded.bubbles.gz vcfs/diploid.shaded
cmd gzip -d -f vcfs/diploid.shaded.vcf.gz

cmd time $THREADCTX -t 1 $SELIST $NUM_INDIVS graphs/pop.se.ctp graphs/pop.ctx graphs/pop.ctx
cmd time $THREADCTX -t 1 $PELIST $NUM_INDIVS graphs/pop.pe.ctp graphs/pop.ctx graphs/pop.ctx
cmd time $THREADCTX -t 1 $SEPELIST $NUM_INDIVS graphs/pop.sepe.ctp graphs/pop.ctx graphs/pop.ctx
cmd mv gap_sizes.* mp_sizes.* graphs/

for x in se pe sepe
do
  cmd time $CALLCTX -t 1 -m 100MB -p graphs/pop.$x.ctp graphs/pop.ctx bubbles/diploid.$x.bubbles.gz
  cmd $PROCCTX bubbles/diploid.$x.bubbles.gz vcfs/diploid.$x
  cmd gzip -d -f vcfs/diploid.$x.vcf.gz
done

# G is a list of genome fasta files
# MG is list of Mask+Genome files
GLIST=
MGLIST=
for i in $(seq 0 $LASTCHROM)
do
  GLIST="$GLIST genomes/genome$i.fa"
  MGLIST="$MGLIST genomes/genome$i.fa genomes/mask$i.fa"
done

# Generate truth VCF
cmd "$BIOINF/sim_mutations/sim_vcf.pl $KMER $MGLIST > vcfs/truth.vcf"

# Compare
msg == Released Cortex ==
LVCF=vcfs/truth.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF vcfs/diploid.oldbc.vcf vcfs/truth.oldbc.vcf OLDBC vcfs/falsepos.oldbc.vcf $GLIST
cmd $HAPLEN vcfs/diploid.oldbc.vcf
msg == New Bubble Caller ==
LVCF=vcfs/truth.oldbc.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF vcfs/diploid.newbc.vcf vcfs/truth.newbc.vcf NEWBC vcfs/falsepos.newbc.vcf $GLIST
cmd $HAPLEN vcfs/diploid.newbc.vcf
msg == Shaded Bubble Caller ==
LVCF=vcfs/truth.newbc.vcf
cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF vcfs/diploid.shaded.vcf vcfs/truth.shaded.vcf SHADED vcfs/falsepos.shaded.vcf $GLIST
cmd $HAPLEN vcfs/diploid.shaded.vcf
LVCF=vcfs/truth.shaded.vcf

for x in se pe sepe
do
  msg == Paths $x ==
  cmd $BIOINF/sim_mutations/sim_compare.pl $LVCF vcfs/diploid.$x.vcf vcfs/truth.$x.vcf PAC vcfs/falsepos.$x.vcf $GLIST
  cmd $HAPLEN vcfs/diploid.$x.vcf
  LVCF=vcfs/truth.$x.vcf
done

msg == Truth ==
cmd $HAPLEN vcfs/truth.vcf

# exit

# Map 5' flanks
STAMPY=stampy.py

if [ $(uname -s) == "Darwin" ]; then
  PYTHON="python2.6"
else
  PYTHON="python"
fi

if [ ! -e stampy.py ]; then
  STAMPY="$PYTHON $HOME/bioinf/stampy-1.0.20/stampy.py"
fi

if [ ! -e $STAMPY_HSH.stidx ]; then cmd $STAMPY -G $STAMPY_HSH $INPUT_SEQ; fi
if [ ! -e $STAMPY_HSH.sthash ]; then cmd $STAMPY -g $STAMPY_HSH -H $STAMPY_HSH; fi

cmd "$STAMPY -g $STAMPY_HSH -h $STAMPY_HSH --inputformat=fasta -M vcfs/diploid.oldbc.5pflanks.fa.gz > vcfs/diploid.oldbc.5pflanks.sam"
cmd "$STAMPY -g $STAMPY_HSH -h $STAMPY_HSH --inputformat=fasta -M vcfs/diploid.newbc.5pflanks.fa.gz > vcfs/diploid.newbc.5pflanks.sam"
cmd "$STAMPY -g $STAMPY_HSH -h $STAMPY_HSH --inputformat=fasta -M vcfs/diploid.shaded.5pflanks.fa.gz > vcfs/diploid.shaded.5pflanks.sam"
cmd "$STAMPY -g $STAMPY_HSH -h $STAMPY_HSH --inputformat=fasta -M vcfs/diploid.se.5pflanks.fa.gz > vcfs/diploid.se.5pflanks.sam"
cmd "$STAMPY -g $STAMPY_HSH -h $STAMPY_HSH --inputformat=fasta -M vcfs/diploid.pe.5pflanks.fa.gz > vcfs/diploid.pe.5pflanks.sam"
cmd "$STAMPY -g $STAMPY_HSH -h $STAMPY_HSH --inputformat=fasta -M vcfs/diploid.sepe.5pflanks.fa.gz > vcfs/diploid.sepe.5pflanks.sam"

# Place calls
cmd "time $PLACECTX vcfs/diploid.oldbc.vcf vcfs/diploid.oldbc.5pflanks.sam $INPUT_SEQ > vcfs/diploid.oldbc.decomp.vcf"
cmd "time $PLACECTX vcfs/diploid.newbc.vcf vcfs/diploid.newbc.5pflanks.sam $INPUT_SEQ > vcfs/diploid.newbc.decomp.vcf"
cmd "time $PLACECTX vcfs/diploid.shaded.vcf vcfs/diploid.shaded.5pflanks.sam $INPUT_SEQ > vcfs/diploid.shaded.decomp.vcf"
cmd "time $PLACECTX vcfs/diploid.se.vcf vcfs/diploid.se.5pflanks.sam $INPUT_SEQ > vcfs/diploid.se.decomp.vcf"
cmd "time $PLACECTX vcfs/diploid.pe.vcf vcfs/diploid.pe.5pflanks.sam $INPUT_SEQ > vcfs/diploid.pe.decomp.vcf"
cmd "time $PLACECTX vcfs/diploid.sepe.vcf vcfs/diploid.sepe.5pflanks.sam $INPUT_SEQ > vcfs/diploid.sepe.decomp.vcf"

# Filter, sort calls
cmd "time vcf-pass vcfs/diploid.oldbc.decomp.vcf | vcf-sort > vcfs/diploid.oldbc.decomp.sort.vcf"
cmd "time vcf-pass vcfs/diploid.newbc.decomp.vcf | vcf-sort > vcfs/diploid.newbc.decomp.sort.vcf"
cmd "time vcf-pass vcfs/diploid.shaded.decomp.vcf | vcf-sort > vcfs/diploid.shaded.decomp.sort.vcf"
cmd "time vcf-pass vcfs/diploid.se.decomp.vcf | vcf-sort > vcfs/diploid.se.decomp.sort.vcf"
cmd "time vcf-pass vcfs/diploid.pe.decomp.vcf | vcf-sort > vcfs/diploid.pe.decomp.sort.vcf"
cmd "time vcf-pass vcfs/diploid.sepe.decomp.vcf | vcf-sort > vcfs/diploid.sepe.decomp.sort.vcf"

# Normalise variants (left align)
for f in oldbc newbc shaded se pe sepe
do
  cmd "$BCFTOOLS norm --remove-duplicate -f $INPUT_SEQ vcfs/diploid.$f.decomp.sort.vcf > vcfs/diploid.$f.norm.vcf"
done

# Generate truth decomp VCF
cmd "zcat -f $INPUT_SEQ | $BIOINF/sim_mutations/sim_decomp_vcf.pl - $GLIST > vcfs/truth.decomp.vcf"
cmd "$BCFTOOLS norm --remove-duplicate -f $INPUT_SEQ vcfs/truth.decomp.vcf > vcfs/truth.norm.vcf"

# Traversal statistics

cmd time $TRAVERSE --colour 0 graphs/pop.ctx

for x in se pe sepe
do
  cmd time $TRAVERSE --colour 0 -p graphs/pop.$x.ctp graphs/pop.ctx
done
