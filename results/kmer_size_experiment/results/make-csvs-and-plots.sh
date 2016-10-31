#!/bin/bash
set -eou pipefail

./make-csv.sh ../perfect_cov/k*/stats.plain.txt > perfect.plain.csv
./make-csv.sh ../perfect_cov/k*/stats.links.txt > perfect.links.csv
./make-csv.sh ../perfect_cov/k*/stats.pe.txt > perfect.pe.csv
./plot-n50-and-errs.R perfect.plain.csv perfect.links.csv perfect.pe.csv "Perfect coverage (100X, 100bp reads)" perfect.cov.pdf

./make-csv.sh ../stoch_cov/k*/stats.plain.txt > stoch.plain.csv
./make-csv.sh ../stoch_cov/k*/stats.links.txt > stoch.links.csv
./make-csv.sh ../stoch_cov/k*/stats.pe.txt > stoch.pe.csv
./plot-n50-and-errs.R stoch.plain.csv stoch.links.csv stoch.pe.csv "Stochastic coverage (100X, 100bp reads)" stoch.cov.pdf

./make-csv.sh ../stocherr_cov/k*/stats.plain.txt > stocherr.plain.csv
./make-csv.sh ../stocherr_cov/k*/stats.links.txt > stocherr.links.csv
./make-csv.sh ../stocherr_cov/k*/stats.pe.txt > stocherr.pe.csv
./plot-n50-and-errs.R stocherr.plain.csv stocherr.links.csv stocherr.pe.csv "Stochastic coverage + Error (100X, 100bp reads, 0.5% err)" stocherr.cov.pdf

./make-cleaning-table.py ../stocherr_cov/k*/graph.k*.dist.txt > cleaning.table.csv
