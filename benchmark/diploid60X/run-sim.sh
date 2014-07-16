#!/bin/bash

# Single diploid 1Mb genome 60X covg
#
# To run:
#   cd dir/this/is/in
#   ./run-sim.sh
#
# To clear up:
#   cd dir/this/is/in
#   ./run-sim.sh clean
#

make -f ../calling-comparison.mk \
  SEQ=../chr21.1Mb.fa NUM_INDIVS=1 PLOIDY=2 KMER=31 \
  SNPS=15000 INDELS=7500 INV=750 INVLEN=100 \
  READLEN=100 MPSIZE=250 ALLELECOVG=30 \
  MEMWIDTH=20 MEMHEIGHT=20 MAPARGS='--substitutionrate=0.01 ' $@

# sites=0.01*sum(1/1)*L=0.01*1M = 10,000
# sites=0.01*sum(1/1+1/2)*L = 15,000
