#!/bin/bash
set -eou pipefail

mkdir -p latest
./make-csv.sh ../perfect_cov/k*/stats.plain.txt > latest/perfect.plain.csv
./make-csv.sh ../perfect_cov/k*/stats.links.txt > latest/perfect.links.csv
./make-csv.sh ../perfect_cov/k*/stats.pe.txt > latest/perfect.pe.csv
./plot-n50-and-errs.R "Perfect cov. (100X, 100bp reads)" latest/perfect.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv latest/perfect.pe.csv

./plot-n50-and-errs.R "Perfect cov. (100X, 100bp reads)" latest/perfect_nope.pdf \
  latest/perfect.plain.csv latest/perfect.links.csv

./make-csv.sh ../stoch_cov/k*/stats.plain.txt > latest/stoch.plain.csv
./make-csv.sh ../stoch_cov/k*/stats.links.txt > latest/stoch.links.csv
./make-csv.sh ../stoch_cov/k*/stats.pe.txt > latest/stoch.pe.csv
./plot-n50-and-errs.R "Stochastic cov. (100X, 100bp reads)" latest/stoch.pdf \
  latest/stoch.plain.csv latest/stoch.links.csv latest/stoch.pe.csv

./make-csv.sh ../stocherr_cov/k*/stats.plain.txt > latest/stocherr.plain.csv
./make-csv.sh ../stocherr_cov/k*/stats.links.txt > latest/stocherr.links.csv
./make-csv.sh ../stocherr_cov/k*/stats.pe.txt > latest/stocherr.pe.csv
./plot-n50-and-errs.R "Stochastic cov. + 0.5% err (100X, 100bp reads)" latest/stocherr.pdf \
  latest/stocherr.plain.csv latest/stocherr.links.csv latest/stocherr.pe.csv

./make-csv.sh ../stocherr_corr/k*/stats.plain.txt > latest/stocherrcorr.plain.csv
./make-csv.sh ../stocherr_corr/k*/stats.links.txt > latest/stocherrcorr.links.csv
./make-csv.sh ../stocherr_corr/k*/stats.pe.txt > latest/stocherrcorr.pe.csv
./plot-bfc.R "Stochastic cov. + 0.5% err (100X, 100bp reads)" latest/stocherrcorr.pdf \
  latest/stocherr.plain.csv latest/stocherr.links.csv latest/stocherr.pe.csv \
  latest/stocherrcorr.plain.csv latest/stocherrcorr.links.csv latest/stocherrcorr.pe.csv

./make-cleaning-table.py ../stocherr_cov/k*/graph.k*.dist.txt > latest/cleaning.table.csv
./make-cleaning-table.py ../stocherr_corr/k*/graph.k*.dist.txt > latest/cleaning.corr.table.csv
