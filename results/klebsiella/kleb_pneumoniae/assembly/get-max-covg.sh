#!/bin/bash

# usage: get-max-covg.sh <genome-size> <in.fa>
# output: <number_of_contigs> <sum_of_lengths> <longest_len> <shortest_len>

genome=$1
seqfile=$2

dnacat -L "$seqfile" | cut -f2 | sort -rn | \
  awk '{if(x+$1>'$genome'){exit;} x+=$1; n+=1; l=$1; if(!f){f=$1}} END{print n,x,f,l}'
