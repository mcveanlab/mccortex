#!/bin/bash

set -eou pipefail
set -x

# Adapted from SGA:
#  https://github.com/jts/sga/blob/726e2e2831796af63fafcf7b52dc80329fa7f31d/src/examples/sga-ecoli-miseq.sh

if [ $# -ne 3 ]
then
  >&2 echo "usage: $0 <in.1.fq.gz> <in.2.fq.gz> <kmers>"
  exit -1
fi

IN1="$1"
IN2="$2"
KMERS="$3"

#
# Parameters
#

# Program paths
CTXDIR=../../
BWA_BIN=$CTXDIR/libs/bwa/bwa
SAMTOOLS_BIN=$CTXDIR/libs/samtools/samtools
SGA_BIN=sga
BAM2DE_BIN=sga-bam2de.pl
ASTAT_BIN=sga-astat.py

# The number of threads to use
CPU=1

# Correction k-mer
CORRECTION_K=41

# Branch trim length
TRIM_LENGTH=100

#
# Dependency checks
#

# Check the required programs are installed and executable
prog_list="$SGA_BIN"
for prog in $prog_list; do
    hash $prog 2>/dev/null || { echo "Error $prog not found. Please place $prog on your PATH or update the *_BIN variables in this script"; exit 1; }
done 

# Check the files are found
file_list="$IN1 $IN2"
for input in $file_list; do
    if [ ! -f $input ]; then
        >&2 echo "Error input file $input not found"; exit 1;
    fi
done

#
# Preprocessing
#

# Preprocess the data to remove ambiguous basecalls
$SGA_BIN preprocess --pe-mode 1 -o reads.pp.fastq $IN1 $IN2

#
# Error Correction
#

# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
$SGA_BIN index -a ropebwt -t $CPU --no-reverse reads.pp.fastq

# Perform k-mer based error correction.
# The k-mer cutoff parameter is learned automatically.
# On Mac only works with 1 thread (-t 1)
$SGA_BIN correct -k $CORRECTION_K --learn -t 1 -o reads.ec.fastq reads.pp.fastq

#
# Primary (contig) assembly
#

# Index the corrected data.
$SGA_BIN index -a ropebwt -t $CPU reads.ec.fastq

# Remove exact-match duplicates and reads with low-frequency k-mers
$SGA_BIN filter -x 2 -t $CPU reads.ec.fastq

for K in `echo "$KMERS"`
do
  mkdir -p k$K
  cd k$K

  # K is the minimum overlap to use when computing the graph.
  # The final assembly can be performed with this overlap or greater

  # Compute the structure of the string graph
  $SGA_BIN overlap -m $K -t $CPU ../reads.ec.filter.pass.fa

  # Perform the contig assembly
  $SGA_BIN assemble -m $K -o assemble.m$K reads.ec.filter.pass.asqg.gz

  cd ..
done
