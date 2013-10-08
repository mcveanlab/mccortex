#!/bin/bash

# Ten diploid 1Kb genomes 20X covg
#
# To run:
#   ./run-sim.sh
#
# To clear up:
#   ./run-sim.sh clean
#

make -f ../calling-comparison.mk \
  SEQ=smaller.fa NUM_INDIVS=10 PLOIDY=2 KMER=31 \
  SNPS=10 INDELS=10 INV=10 INVLEN=10 \
  READLEN=100 MPSIZE=250 ALLELECOVG=10 \
  MEMWIDTH=20 MEMHEIGHT=15 $@
