./make-csv.sh ../perfect_cov/k*/stats.plain.txt > perfect.plain.csv
./make-csv.sh ../perfect_cov/k*/stats.links.txt > perfect.links.csv
./plot-n50-and-errs.R perfect.plain.csv perfect.links.csv "Perfect coverage (100X, 100bp reads)" perfect.cov.pdf

./make-csv.sh ../stoch_cov/k*/stats.plain.txt > stoch.plain.csv
./make-csv.sh ../stoch_cov/k*/stats.links.txt > stoch.links.csv
./plot-n50-and-errs.R stoch.plain.csv stoch.links.csv "Stochastic coverage (100X, 100bp reads)" stoch.cov.pdf

./make-csv.sh ../stocherr_cov/k*/stats.plain.txt > stocherr.plain.csv
./make-csv.sh ../stocherr_cov/k*/stats.links.txt > stocherr.links.csv
./plot-n50-and-errs.R stocherr.plain.csv stocherr.links.csv "Stochastic coverage + Error (100X, 100bp reads, 0.5% err)" stocherr.cov.pdf
