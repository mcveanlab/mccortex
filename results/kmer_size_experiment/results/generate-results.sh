#!/bin/bash
set -eou pipefail

# hide k99 files
echo "-- hiding k=99"
for d in perfect_cov stoch_cov stocherr_cov stocherr_corr
do
  if [ -d ../$d/k99 ]
  then
    mv ../$d/k99 ../$d/hidden_k99
  fi
done

mkdir -p latest
echo "-- Perfect Coverage"
./make-csv.sh ../perfect_cov/k*/stats.plain.txt > latest/perfect.plain.csv
./make-csv.sh ../perfect_cov/k*/stats.links.txt > latest/perfect.links.csv
./make-csv.sh ../perfect_cov/k*/stats.pe.txt > latest/perfect.pe.csv
./plot-ng50-and-errs.R "Perfect cov. (100X, 100bp reads)" latest/perfect.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv latest/perfect.pe.csv

./plot-ng50-and-errs.R "Perfect cov. (100X, 100bp reads)" latest/perfect_no_pe.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv

echo "-- Stochastic Coverage"
./make-csv.sh ../stoch_cov/k*/stats.plain.txt > latest/stoch.plain.csv
./make-csv.sh ../stoch_cov/k*/stats.links.txt > latest/stoch.links.csv
./make-csv.sh ../stoch_cov/k*/stats.pe.txt > latest/stoch.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. (100X, 100bp reads)" latest/stoch.pdf \
  latest/stoch.plain.csv latest/stoch.links.csv latest/stoch.pe.csv

echo "-- Stochastic Coverage + Error"
./make-csv.sh ../stocherr_cov/k*/stats.plain.txt > latest/stocherr.plain.csv
./make-csv.sh ../stocherr_cov/k*/stats.links.txt > latest/stocherr.links.csv
./make-csv.sh ../stocherr_cov/k*/stats.pe.txt > latest/stocherr.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. + 0.5% err (100X, 100bp reads)" latest/stocherr.pdf \
  latest/stocherr.plain.csv latest/stocherr.links.csv latest/stocherr.pe.csv

echo "-- Stochastic Coverage + Error + Error Correction"
./make-csv.sh ../stocherr_corr/k*/stats.plain.txt > latest/stocherrcorr.plain.csv
./make-csv.sh ../stocherr_corr/k*/stats.links.txt > latest/stocherrcorr.links.csv
./make-csv.sh ../stocherr_corr/k*/stats.pe.txt > latest/stocherrcorr.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. + 0.5% err + correct (100X, 100bp reads)" latest/stocherrcorr.pdf \
  latest/stocherr.plain.csv latest/stocherr.links.csv latest/stocherr.pe.csv \
  latest/stocherrcorr.plain.csv latest/stocherrcorr.links.csv latest/stocherrcorr.pe.csv

# Gather SGA results
if [ -d ../stocherr_cov/sga ]
then
  echo "-- SGA plots"
  # ../stocherr_cov/sga/k21/stats.k21.txt
  ./make-csv.sh ../stocherr_cov/sga/k*/stats.k*.txt > latest/stocherr.sga.csv
  ./plot-mccortex-vs-sga.R latest/links-vs-sga-ng50.pdf latest/links-vs-sga-errs.pdf latest/stocherr.links.csv latest/stocherr.sga.csv
  ./plot-mccortex-vs-sga.R latest/pe-vs-sga-ng50.pdf latest/pe-vs-sga-errs.pdf latest/stocherr.pe.csv latest/stocherr.sga.csv
fi

if [ -d ../stocherr_corr/sga ]
then
  echo "-- SGA + Error Correction plots"
  # ../stocherr_corr/sga/k21/stats.k21.txt
  ./make-csv.sh ../stocherr_corr/sga/k*/stats.k*.txt > latest/stocherrcorr.sga.csv
  ./plot-mccortex-vs-sga.R latest/corr-links-vs-sga-ng50.pdf latest/corr-links-vs-sga-errs.pdf latest/stocherrcorr.links.csv latest/stocherrcorr.sga.csv
  ./plot-mccortex-vs-sga.R latest/corr-pe-vs-sga-ng50.pdf latest/corr-pe-vs-sga-errs.pdf latest/stocherrcorr.pe.csv latest/stocherrcorr.sga.csv
fi

echo "-- Plain vs links"
./plot-ng50-three-sets.R latest/plain-vs-links.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv \
  latest/stoch.plain.csv latest/stoch.links.csv \
  latest/stocherr.plain.csv latest/stocherr.links.csv

echo "-- Plain vs PE"
./plot-ng50-three-sets.R latest/plain-vs-pe.pdf \
  latest/perfect.plain.csv latest/perfect.pe.csv \
  latest/stoch.plain.csv latest/stoch.pe.csv \
  latest/stocherr.plain.csv latest/stocherr.pe.csv

echo "-- Plain vs links (with error correction)"
./plot-ng50-three-sets.R latest/plain-vs-links-corr.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv \
  latest/stoch.plain.csv latest/stoch.links.csv \
  latest/stocherrcorr.plain.csv latest/stocherrcorr.links.csv

echo "-- Plain vs PE links (with error correction)"
./plot-ng50-three-sets.R latest/plain-vs-pe-corr.pdf \
  latest/perfect.plain.csv latest/perfect.pe.csv \
  latest/stoch.plain.csv latest/stoch.pe.csv \
  latest/stocherrcorr.plain.csv latest/stocherrcorr.pe.csv

echo "-- Making cleaning tables"
./make-cleaning-table.py ../stocherr_cov/k*/graph.k*.dist.txt > latest/cleaning.table.csv
./make-cleaning-table.py ../stocherr_corr/k*/graph.k*.dist.txt > latest/cleaning.corr.table.csv

echo "-- Make link count csv"
for t in se pe; do
  cat ../perfect_cov/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > latest/perfect.linkcounts.$t.csv
  cat ../stoch_cov/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > latest/stoch.linkcounts.$t.csv
  cat ../stocherr_cov/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > latest/stocherr.linkcounts.$t.csv
  cat ../stocherr_corr/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > latest/stocherrcorr.linkcounts.$t.csv
done
for t in se pe; do
  for s in perfect stoch stocherr stocherrcorr; do
    ./plot-link-counts.R latest/$s.linkcounts.$t.pdf latest/$s.linkcounts.$t.csv
  done
  ./plot-link-counts-together.R latest/linkcounts.$t.pdf latest/{perfect,stoch,stocherr,stocherrcorr}.linkcounts.$t.csv
done


# unhide k99 files
echo "-- recovering k=99"
for d in perfect_cov stoch_cov stocherr_cov stocherr_corr
do
  if [ -d ../$d/hidden_k99 ]
  then
    mv ../$d/hidden_k99 ../$d/k99
  fi
done
