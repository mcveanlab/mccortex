#!/bin/bash
set -eou pipefail
set -x

# Set to 1 to filter out k=99, else include k99
HIDE_K99=0

if (($HIDE_K99)); then
  OUTDIR=latest_k91
else
  OUTDIR=latest_k99
fi

mkdir -p $OUTDIR

# hide k99 files
if (($HIDE_K99)); then
  echo "-- hiding k=99"
  for d in perfect_cov stoch_cov stocherr_cov stocherr_corr
  do
    if [ -d ../$d/k99 ]
    then
      mv ../$d/k99 ../$d/hidden_k99
    fi
  done
fi

echo "-- Perfect Coverage"
./make-csv.sh ../perfect_cov/k*/stats.plain.txt > $OUTDIR/perfect.plain.csv
./make-csv.sh ../perfect_cov/k*/stats.links.txt > $OUTDIR/perfect.links.csv
./make-csv.sh ../perfect_cov/k*/stats.pe.txt > $OUTDIR/perfect.pe.csv
./plot-ng50-and-errs.R "Perfect cov. (100X, 100bp reads)" $OUTDIR/perfect.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.links.csv $OUTDIR/perfect.pe.csv

./plot-ng50-and-errs.R "Perfect cov. (100X, 100bp reads)" $OUTDIR/perfect_no_pe.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.links.csv

echo "-- Stochastic Coverage"
./make-csv.sh ../stoch_cov/k*/stats.plain.txt > $OUTDIR/stoch.plain.csv
./make-csv.sh ../stoch_cov/k*/stats.links.txt > $OUTDIR/stoch.links.csv
./make-csv.sh ../stoch_cov/k*/stats.pe.txt > $OUTDIR/stoch.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. (100X, 100bp reads)" $OUTDIR/stoch.pdf \
  $OUTDIR/stoch.plain.csv $OUTDIR/stoch.links.csv $OUTDIR/stoch.pe.csv

echo "-- Stochastic Coverage + Error"
./make-csv.sh ../stocherr_cov/k*/stats.plain.txt > $OUTDIR/stocherr.plain.csv
./make-csv.sh ../stocherr_cov/k*/stats.links.txt > $OUTDIR/stocherr.links.csv
./make-csv.sh ../stocherr_cov/k*/stats.pe.txt > $OUTDIR/stocherr.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. + 0.5% err (100X, 100bp reads)" $OUTDIR/stocherr.pdf \
  $OUTDIR/stocherr.plain.csv $OUTDIR/stocherr.links.csv $OUTDIR/stocherr.pe.csv

echo "-- Stochastic Coverage + Error + Error Correction"
./make-csv.sh ../stocherr_corr/k*/stats.plain.txt > $OUTDIR/stocherrcorr.plain.csv
./make-csv.sh ../stocherr_corr/k*/stats.links.txt > $OUTDIR/stocherrcorr.links.csv
./make-csv.sh ../stocherr_corr/k*/stats.pe.txt > $OUTDIR/stocherrcorr.pe.csv
./plot-ng50-and-errs.R "Stochastic cov. + 0.5% err + correct (100X, 100bp reads)" $OUTDIR/stocherrcorr.pdf \
  $OUTDIR/stocherr.plain.csv $OUTDIR/stocherr.links.csv $OUTDIR/stocherr.pe.csv \
  $OUTDIR/stocherrcorr.plain.csv $OUTDIR/stocherrcorr.links.csv $OUTDIR/stocherrcorr.pe.csv

# Gather SGA results
if [ -d ../stocherr_cov/sga ]
then
  echo "-- SGA plots"
  # ../stocherr_cov/sga/k21/stats.k21.txt
  ./make-csv.sh ../stocherr_cov/sga/k*/stats.k*.txt > $OUTDIR/stocherr.sga.csv
  ./plot-mccortex-vs-sga.R $OUTDIR/links-vs-sga-ng50.pdf $OUTDIR/links-vs-sga-errs.pdf $OUTDIR/stocherr.links.csv $OUTDIR/stocherr.sga.csv
  ./plot-mccortex-vs-sga.R $OUTDIR/pe-vs-sga-ng50.pdf $OUTDIR/pe-vs-sga-errs.pdf $OUTDIR/stocherr.pe.csv $OUTDIR/stocherr.sga.csv
fi

if [ -d ../stocherr_corr/sga ]
then
  echo "-- SGA + Error Correction plots"
  # ../stocherr_corr/sga/k21/stats.k21.txt
  ./make-csv.sh ../stocherr_corr/sga/k*/stats.k*.txt > $OUTDIR/stocherrcorr.sga.csv
  ./plot-mccortex-vs-sga.R $OUTDIR/corr-links-vs-sga-ng50.pdf $OUTDIR/corr-links-vs-sga-errs.pdf $OUTDIR/stocherrcorr.links.csv $OUTDIR/stocherrcorr.sga.csv
  ./plot-mccortex-vs-sga.R $OUTDIR/corr-pe-vs-sga-ng50.pdf $OUTDIR/corr-pe-vs-sga-errs.pdf $OUTDIR/stocherrcorr.pe.csv $OUTDIR/stocherrcorr.sga.csv
  # McCortex (SE) vs McCortex (PE) vs SGA (all with bfc corrected reads)
  ./plot-mccortex-se-pe-vs-sga.R $OUTDIR/corr-both-vs-sga-ng50.pdf $OUTDIR/corr-both-vs-sga-errs.pdf $OUTDIR/stocherrcorr.links.csv $OUTDIR/stocherrcorr.pe.csv $OUTDIR/stocherrcorr.sga.csv
  # McCortex (SE) vs McCortex (PE) vs SGA (not with bfc corrected reads)
  ./plot-mccortex-se-pe-vs-sga.R $OUTDIR/both-vs-sga-ng50.pdf $OUTDIR/both-vs-sga-errs.pdf $OUTDIR/stocherr.links.csv $OUTDIR/stocherr.pe.csv $OUTDIR/stocherr.sga.csv
  # Corrected + McCortex vs raw + SGA
  ./plot-mccortex-vs-sga.R $OUTDIR/corr-links-vs-raw-sga-ng50.pdf $OUTDIR/corr-links-vs-raw-sga-errs.pdf $OUTDIR/stocherrcorr.links.csv $OUTDIR/stocherr.sga.csv
  ./plot-mccortex-vs-sga.R $OUTDIR/corr-pe-vs-raw-sga-ng50.pdf $OUTDIR/corr-pe-vs-raw-sga-errs.pdf $OUTDIR/stocherrcorr.pe.csv $OUTDIR/stocherr.sga.csv
fi

echo "-- Plain vs links"
./plot-ng50-three-sets.R $OUTDIR/plain-vs-links-ng50.pdf $OUTDIR/plain-vs-links-errs.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.links.csv \
  $OUTDIR/stoch.plain.csv $OUTDIR/stoch.links.csv \
  $OUTDIR/stocherr.plain.csv $OUTDIR/stocherr.links.csv

echo "-- Plain vs PE"
./plot-ng50-three-sets.R $OUTDIR/plain-vs-pe-ng50.pdf $OUTDIR/plain-vs-pe-errs.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.pe.csv \
  $OUTDIR/stoch.plain.csv $OUTDIR/stoch.pe.csv \
  $OUTDIR/stocherr.plain.csv $OUTDIR/stocherr.pe.csv

echo "-- Plain vs links (with error correction)"
./plot-ng50-three-sets.R $OUTDIR/plain-vs-links-corr-ng50.pdf $OUTDIR/plain-vs-links-corr-errs.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.links.csv \
  $OUTDIR/stoch.plain.csv $OUTDIR/stoch.links.csv \
  $OUTDIR/stocherrcorr.plain.csv $OUTDIR/stocherrcorr.links.csv

echo "-- Plain vs PE links (with error correction)"
./plot-ng50-three-sets.R $OUTDIR/plain-vs-pe-corr-ng50.pdf $OUTDIR/plain-vs-pe-corr-errs.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.pe.csv \
  $OUTDIR/stoch.plain.csv $OUTDIR/stoch.pe.csv \
  $OUTDIR/stocherrcorr.plain.csv $OUTDIR/stocherrcorr.pe.csv

echo "-- Plain vs links (perfect, stochastic, error, corrected)"
./plot-ng50-four-sets.R $OUTDIR/plain-vs-links-fourway-ng50.pdf $OUTDIR/plain-vs-links-fourway-errs.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.links.csv \
  $OUTDIR/stoch.plain.csv $OUTDIR/stoch.links.csv \
  $OUTDIR/stocherr.plain.csv $OUTDIR/stocherr.links.csv \
  $OUTDIR/stocherrcorr.plain.csv $OUTDIR/stocherrcorr.links.csv

echo "-- Plain vs PE (perfect, stochastic, error, corrected)"
./plot-ng50-four-sets.R $OUTDIR/plain-vs-pe-fourway-ng50.pdf $OUTDIR/plain-vs-pe-fourway-errs.pdf \
  $OUTDIR/perfect.plain.csv $OUTDIR/perfect.pe.csv \
  $OUTDIR/stoch.plain.csv $OUTDIR/stoch.pe.csv \
  $OUTDIR/stocherr.plain.csv $OUTDIR/stocherr.pe.csv \
  $OUTDIR/stocherrcorr.plain.csv $OUTDIR/stocherrcorr.pe.csv

echo "-- Making cleaning tables"
./make-cleaning-table.py ../stocherr_cov/k*/graph.k*.dist.txt > $OUTDIR/cleaning.table.csv
./make-cleaning-table.py ../stocherr_corr/k*/graph.k*.dist.txt > $OUTDIR/cleaning.corr.table.csv

echo "-- Make link count csv"
for t in se pe; do
  cat ../perfect_cov/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > $OUTDIR/perfect.linkcounts.$t.csv
  cat ../stoch_cov/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > $OUTDIR/stoch.linkcounts.$t.csv
  cat ../stocherr_cov/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > $OUTDIR/stocherr.linkcounts.$t.csv
  cat ../stocherr_corr/k*/graph.k*.$t.raw.ctp.gz.log | ./count-links.pl > $OUTDIR/stocherrcorr.linkcounts.$t.csv
done
for t in se pe; do
  for s in perfect stoch stocherr stocherrcorr; do
    ./plot-link-counts.R "$OUTDIR/$s"_linkcounts_$t.pdf $OUTDIR/$s.linkcounts.$t.csv
  done
  ./plot-link-counts-threeway.R $OUTDIR/linkcounts_$t\_threeway.pdf $OUTDIR/{perfect,stoch,stocherr}.linkcounts.$t.csv
  ./plot-link-counts-together.R $OUTDIR/linkcounts_$t\_fourway.pdf $OUTDIR/{perfect,stoch,stocherr,stocherrcorr}.linkcounts.$t.csv
done


# unhide k99 files
if (($HIDE_K99)); then
  echo "-- recovering k=99"
  for d in perfect_cov stoch_cov stocherr_cov stocherr_corr
  do
    if [ -d ../$d/hidden_k99 ]
    then
      mv ../$d/hidden_k99 ../$d/k99
    fi
  done
fi
