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
  SNPS=10000 INDELS=5000 INV=500 INVLEN=100 \
  READLEN=100 MPSIZE=250 ALLELECOVG=30 \
  MEMWIDTH=20 MEMHEIGHT=20 $@
