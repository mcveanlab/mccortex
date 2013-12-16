#!/bin/bash

# Single small diploid 1Kb 20X covg
#
# To run:
#   ./run-sim.sh
#
# To clear up:
#   ./run-sim.sh clean
#

{ echo '>rnd'; seqrnd 1000; } | facat > rnd.fa

make -f ../calling-comparison.mk \
  SEQ=rnd.fa NUM_INDIVS=1 PLOIDY=2 KMER=31 \
  SNPS=0 INDELS=100 INV=0 INVLEN=10 \
  READLEN=100 MPSIZE=250 ALLELECOVG=10 \
  MEMWIDTH=20 MEMHEIGHT=15 $@
