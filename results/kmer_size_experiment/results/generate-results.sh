#!/bin/bash
set -eou pipefail

# hide k99 files
for d in perfect_cov stoch_cov stocherr_cov stocherr_corr
do
  if [ -d ../$d/k99 ]
  then
    mv ../$d/k99 ../$d/99_hidden
  fi
done

mkdir -p latest
./make-csv.sh ../perfect_cov/k*/stats.plain.txt > latest/perfect.plain.csv
./make-csv.sh ../perfect_cov/k*/stats.links.txt > latest/perfect.links.csv
./make-csv.sh ../perfect_cov/k*/stats.pe.txt > latest/perfect.pe.csv
./plot-ng50-and-errs.R "Perfect cov. (100X, 100bp reads)" latest/perfect.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv latest/perfect.pe.csv

./plot-ng50-and-errs.R "Perfect cov. (100X, 100bp reads)" latest/perfect_no_pe.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv

./make-csv.sh ../stoch_cov/k*/stats.plain.txt > latest/stoch.plain.csv
./make-csv.sh ../stoch_cov/k*/stats.links.txt > latest/stoch.links.csv
./make-csv.sh ../stoch_cov/k*/stats.pe.txt > latest/stoch.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. (100X, 100bp reads)" latest/stoch.pdf \
  latest/stoch.plain.csv latest/stoch.links.csv latest/stoch.pe.csv

./make-csv.sh ../stocherr_cov/k*/stats.plain.txt > latest/stocherr.plain.csv
./make-csv.sh ../stocherr_cov/k*/stats.links.txt > latest/stocherr.links.csv
./make-csv.sh ../stocherr_cov/k*/stats.pe.txt > latest/stocherr.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. + 0.5% err (100X, 100bp reads)" latest/stocherr.pdf \
  latest/stocherr.plain.csv latest/stocherr.links.csv latest/stocherr.pe.csv

./make-csv.sh ../stocherr_corr/k*/stats.plain.txt > latest/stocherrcorr.plain.csv
./make-csv.sh ../stocherr_corr/k*/stats.links.txt > latest/stocherrcorr.links.csv
./make-csv.sh ../stocherr_corr/k*/stats.pe.txt > latest/stocherrcorr.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. + 0.5% err + correct (100X, 100bp reads)" latest/stocherrcorr.pdf \
  latest/stocherr.plain.csv latest/stocherr.links.csv latest/stocherr.pe.csv \
  latest/stocherrcorr.plain.csv latest/stocherrcorr.links.csv latest/stocherrcorr.pe.csv

# Gather SGA results
if [ -d ../stocherr_cov/sga ]
then
  # ../stocherr_cov/sga/k21/stats.k21.txt
  ./make-csv.sh ../stocherr_cov/sga/k*/stats.k*.txt > latest/stocherr.sga.csv
  ./plot-mccortex-vs-sga.R latest/links-vs-sga-ng50.pdf latest/links-vs-sga-errs.pdf latest/stocherr.links.csv latest/stocherr.sga.csv
  ./plot-mccortex-vs-sga.R latest/pe-vs-sga-ng50.pdf latest/pe-vs-sga-errs.pdf latest/stocherr.pe.csv latest/stocherr.sga.csv
fi

./plot-ng50-three-sets.R latest/plain-vs-links.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv \
  latest/stoch.plain.csv latest/stoch.links.csv \
  latest/stocherr.plain.csv latest/stocherr.links.csv

./plot-ng50-three-sets.R latest/plain-vs-pe.pdf \
  latest/perfect.plain.csv latest/perfect.pe.csv \
  latest/stoch.plain.csv latest/stoch.pe.csv \
  latest/stocherr.plain.csv latest/stocherr.pe.csv

./plot-ng50-three-sets.R latest/plain-vs-links-corr.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv \
  latest/stoch.plain.csv latest/stoch.links.csv \
  latest/stocherrcorr.plain.csv latest/stocherrcorr.links.csv

./plot-ng50-three-sets.R latest/plain-vs-pe-corr.pdf \
  latest/perfect.plain.csv latest/perfect.pe.csv \
  latest/stoch.plain.csv latest/stoch.pe.csv \
  latest/stocherrcorr.plain.csv latest/stocherrcorr.pe.csv

./make-cleaning-table.py ../stocherr_cov/k*/graph.k*.dist.txt > latest/cleaning.table.csv
./make-cleaning-table.py ../stocherr_corr/k*/graph.k*.dist.txt > latest/cleaning.corr.table.csv

# unhide k99 files
for d in perfect_cov stoch_cov stocherr_cov stocherr_corr
do
  if [ -d ../$d/k99_hidden ]
  then
    mv ../$d/99_hidden ../$d/k99
  fi
done
