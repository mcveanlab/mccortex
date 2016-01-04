#!/bin/bash

set -eou pipefail

snpfile=$1
reffile=$2

echo '##fileformat=VCFv4.1'
echo '##fileDate='`date '+%Y%m%d'`
echo "##reference=$reffile"
~/c/dnacat/bin/dnacat -L $reffile | awk '{OFS=""; print "##chrom=<ID=",$1,",length=",$2,">"}'
echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
echo | awk '{OFS="\t"; print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"}'

grep -E '^\s*[0-9]' $snpfile | awk '{OFS="\t"; print $14,$1,".",$2,$3,".",".",".","GT"}'
